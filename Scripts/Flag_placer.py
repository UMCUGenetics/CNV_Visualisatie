import pysam
from datetime import datetime, date
from statistics import mean, median
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=True, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--threshold', '-t', required=False, default=0, type=int, help='int specifying the minimun amount'
                                                                                    'of reads before a flag is made')
parser.add_argument('--minpercentage', '-mp', required=False, default=0, type=float, help='float specifying the'
                                                                                          'threshold for the minimum'
                                                                                          'percentage of total reads in'
                                                                                          'region before flagged.')
parser.add_argument('--region', '-r', required=False, default='', type=str, help='String specifying the region in '
                                                                                 'format: "chr#:start-stop". use '
                                                                                 'chr#:0-0 for whole chromosome.')
parser.add_argument('--high_insert_size', '-hi', required=False, default=500, type=int,
                    help='Length of insert size to be classified as high.')

args = parser.parse_args()


def fetch_reads():
    """ The fetch_reads function fetches the reads from the bam file

    :return reads: Pysam object containing read information
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')

    chromosome = args.region.split(':')[0].replace('chr', '')
    start = args.region.split(':')[1].split('-')[0]
    end = args.region.split(':')[1].split('-')[1]

    if start == '0' and end == '0':
        reads = bamfile.fetch(chromosome)

    else:
        reads = bamfile.fetch(chromosome, int(start), int(end))

    return reads


def place_flags(reads):
    """ The place_flags function gets the coordinates of interesting regions and returns the coordinates of the flags
    with their information.

    :param reads: Pysam object containing read information
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    all_flags = []
    chromosome = args.region.split(':')[0].replace('chr', '')

    isbuildingflags = [False, False, False, False]
    flags = [[chromosome, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}],
             [chromosome, None, None, {'type': 'high_insert_size', 'count': 0, 'total': 0, 'lengths': []}],
             [chromosome, None, None, {'type': 'unmapped_mate', 'count': 0, 'total': 0}],
             [chromosome, None, None, {'type': 'softclips', 'count': 0, 'total': 0, 'bases': []}]]

    for read in reads:
        if not read.is_unmapped:
            start, end = true_position(read)
            # place flag if read and its mate have the same orientation
            flags, isbuildingflags, all_flags = flag_sameorientation(read, flags, isbuildingflags, all_flags,
                                                                     chromosome, start)

            # place flag if read has high insert size.
            flags, isbuildingflags, all_flags = flag_high_isize(read, flags, isbuildingflags, all_flags, chromosome,
                                                                start)
            # place flag if read has unmapped mate.
            flags, isbuildingflags, all_flags = flag_unmapped_mate(read, flags, isbuildingflags, all_flags, chromosome,
                                                                   start)

            #place flag if read has sofclip bases.
            flags, isbuildingflags, all_flags = flag_softclips(read, flags, isbuildingflags, all_flags, chromosome,
                                                                   start)
            flags = update_total(flags)

    return all_flags


def flag_sameorientation(read, flags, isbuildingflags, all_flags, chromosome, start):
    if issameorientation(read):
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 0)

    elif not issameorientation(read) and isbuildingflags[0] and start > flags[0][2]:
        percentage = round(flags[0][3]['count'] / flags[0][3]['total'], 2)
        if flags[0][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[0])
        flags[0] = [chromosome, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}]
        isbuildingflags[0] = False

    return flags, isbuildingflags, all_flags


def flag_high_isize(read, flags, isbuildingflags, all_flags, chromosome, start):
    if read.isize > args.high_insert_size:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 1)
        flags[1][3]['lengths'].append(read.isize)

    elif not read.isize > args.high_insert_size and isbuildingflags[1] and start > flags[1][2]:
        percentage = round(flags[1][3]['count'] / flags[1][3]['total'], 2)
        if flags[1][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[1])
        flags[1] = [chromosome, None, None, {'type': 'high_insert_size', 'count': 0, 'total': 0, 'lengths': []}]
        isbuildingflags[1] = False

    return flags, isbuildingflags, all_flags


def flag_unmapped_mate(read, flags, isbuildingflags, all_flags, chromosome, start):
    if read.mate_is_unmapped:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 2)

    elif not read.mate_is_unmapped and isbuildingflags[2] and start > flags[2][2]:
        percentage = round(flags[2][3]['count'] / flags[2][3]['total'], 2)
        if flags[2][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[2])
        flags[2] = [chromosome, None, None, {'type': 'unmapped_mate', 'count': 0, 'total': 0}]
        isbuildingflags[2] = False

    return flags, isbuildingflags, all_flags


def flag_softclips(read, flags, isbuildingflags, all_flags, chromosome, start):
    if has_softclips(read):
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 2)
        flags[3][3]['bases'].append(softclipbases(read))

    elif not read.isize > args.high_insert_size and isbuildingflags[3] and start > flags[3][2]:
        percentage = round(flags[3][3]['count'] / flags[3][3]['total'], 3)
        if flags[3][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[3])
        flags[3] = [chromosome, None, None, {'type': 'softclips', 'count': 0, 'total': 0, 'bases': []}]
        isbuildingflags[3] = False

    return flags, isbuildingflags, all_flags


def update_total(flags):
    """ The update_total function iterates over all the flag types and increments the total number of reads it has
    encountered by 1.

    :param flags: a 2d list containing all the flag information.
    return flags: a 2d list containing all the flag information.
    """
    for flag in flags:
        flag[3]['total'] += 1

    return flags


def generate_flag(read, flags, isbuildingflags, flagindex):
    """ The generate_flag function receives a read and decides if it should be included in the current working flag or
    not. Or it starts the creation of a new flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags; a list identifying which types of flags are being built/edited
    :param flagindex: an integer identifying which flag is being built/edited
    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list identifying which types of flags are being built/edited
    """
    start, end = true_position(read)

    if not isbuildingflags[flagindex]:
        flags[flagindex][1] = start
        flags[flagindex][2] = end
        isbuildingflags[flagindex] = True

    if isbuildingflags[flagindex]:
        flags[flagindex][2] = end

    flags[flagindex][3]['count'] += 1

    return flags, isbuildingflags


def true_position(read):
    """ The true_position function receives a read and determines the start of the read as presented in igv by including
    unmapped basepairs.

    :param read: pysam object containing data of a read.
    :return start: an integer indicating the true start of a read.
    :return end: an integer idicating the true end of a read.
    """
    cigar = read.cigar
    start = read.positions[0]
    end = read.positions[-1]

    if cigar[0][0] != 0:
        start -= cigar[0][1]

    if cigar[-1][0] != 0:
        end += cigar[-1][1]

    return start, end


def issameorientation(read):
    """ The issameorientation function returns a bool returning true if a pair of reads have the same orientation and
    false if they have an oposite orientation.

    :param read: Pysam object containing data of a read.
    :return bool: A boolean returning True if the reads of a pair have the same orientation.
    """
    if read.is_paired and not read.mate_is_unmapped:
        if read.is_reverse == read.mate_is_reverse:
            return True

        else:
            return False


def has_softclips(read):
    cigar = read.cigar
    for element in cigar:
        if element[0] == 4:
            return True

    return False


def softclipbases(read):
    cigar = read.cigar
    sofclipbases = 0

    for element in cigar:
        if element[0] == 4:
            sofclipbases += element[1]

    return sofclipbases


def write_bedfile(flags):
    """ The write_bedfile function writes a file in BED format that can be loaded in igv and visualises the read data.

    :param regions: a list of coordinates specified in the vcf file.
    :param read_data: a 2d list containing the read data of every region.
    """
    text = 'track name=same_direction_reads description="Region_Summary." db=hg19 gffTags=on itemRGB="On"\n'

    for flag in flags:
        percentage = round(flag[3]['count'] / flag[3]['total'], 2)

        region = f"{flag[0]}\t{flag[1]}\t{flag[2]}"
        description = f"Name={flag[3]['type']};Readcount={flag[3]['count']};Total_reads={flag[3]['total']};" \
                      f"Percentage={percentage}"

        if flag[3]['type'] == 'high_insert_size':
            lengths = flag[3]['lengths']
            description += f";Avg_insert_size={round(mean(lengths))};Med_insert_size={median(lengths)};" \
                           f"Lower_limit={min(lengths)};Upper_limit={max(lengths)}"

        if flag[3]['type'] == 'softclips':
            bases = flag[3]['bases']
            description += f"Avg_sofclip_bases={round(mean(bases))};Med_sofclip_bases={median(bases)};" \
                           f"Lower_limit={min(bases)};Upper_limit={max(bases)}"

        if flag[3]['type'] == 'same_orientation':
            rgb = f"150,200,150"

        elif flag[3]['type'] == 'high_insert_size':
            rgb = f"200,150,150"

        elif flag[3]['type'] == 'unmapped_mate':
            rgb = f"150,150,200"

        elif flag[3]['type'] == 'sofclips':
            rgb = f"150,200,200"

        text += f"{region}\t{description}\t0\t.\t{flag[1]}\t{flag[2]}\t{rgb}\n"

    with open(args.output + f'/{args.region}_flags.bed', 'w') as bedfile:
        bedfile.write(text)


def write_logfile():
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}\nScript finished at: {current_time} {current_day}\n{"-"*40}' \
           f'Parameters{"-"*40}\nRegion: {args.region}\nBamfile: {args.bam}\nOutput_folder: {args.output}\n' \
           f'Read threshold: {args.threshold}\nInsert size threshold: {args.high_insert_size}\n' \
           f'Minimal percentage: {args.minpercentage}'

    with open(args.output + f'/{args.region}_flag_log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    """ The main calls the other functions in the right order."""

    reads = fetch_reads()  # get all reads at the regions from the vcf call.

    flags = place_flags(reads)  # retrieve data from reads.

    write_bedfile(flags)  # write the result in a BED file.

    if args.log:
        write_logfile()  # write logfile with parameters
