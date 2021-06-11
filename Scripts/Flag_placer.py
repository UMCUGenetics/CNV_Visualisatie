import pysam
from datetime import datetime, date
from statistics import mean, median
import os
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--region', '-r', required=True, type=str,
                    help='String specifying the region in format: "chr#:start-stop". use chr# for whole chromosome.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=False, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--threshold', '-t', required=False, default=0, type=int,
                    help='int specifying the minimun amount of reads before a flag is made')
parser.add_argument('--minpercentage', '-mp', required=False, default=0, type=float,
                    help='float specifying thethreshold for the minimum percentage of total reads in region before '
                         'flagged.')
parser.add_argument('--high_insert_size', '-hi', required=False, default=-1, type=int,
                    help='Length of insert size to be classified as high.')
parser.add_argument('--ultra_high_insert_size', '-uhi', required=False, default=20000, type=int,
                    help='Length of insert size to be classified as ultra high.')
parser.add_argument('--name', '-n', required=False, default='output', type=str,
                    help='Name for the project. This is the name of the output file.')

args = parser.parse_args()


def fetch_reads():
    """ The fetch_reads function fetches the reads from the bam file

    :return reads: Pysam object containing read information
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')

    if ':' in args.region and '-' in args.region:
        chromosome, start, end = re.split(':|-', args.region)
        chromosome = chromosome.replace('chr', '')

        reads = bamfile.fetch(chromosome, int(start), int(end))

    else:
        chromosome = args.region.replace('chr', '')
        reads = bamfile.fetch(chromosome)

    return reads


def place_flags(reads):
    """ The place_flags function gets the coordinates of interesting regions and returns the coordinates of the flags
    with their information.

    :param reads: Pysam object containing read information
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    :return read_data: A list with the number of total reads, unmapped reads and reads with 0 mapped positions.
    """
    all_flags = []
    read_data = [0, 0, 0]  # [0] is number of reads. [1] is unmapped reads. [2] is reads with 0 mapped positions

    isbuildingflags = [False, False, False, False, False, False]
    flags = [[None, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}],
             [None, None, None, {'type': 'high_insert_size', 'count': 0, 'total': 0, 'lengths': []}],
             [None, None, None, {'type': 'unmapped_mate', 'count': 0, 'total': 0}],
             [None, None, None, {'type': 'ultra_high_insert_size', 'count': 0, 'total': 0, 'lengths': []}],
             [None, None, None, {'type': 'inter_chromosomal_pair', 'count': 0, 'total': 0}],
             [None, None, None, {'type': 'face_away', 'count': 0, 'total': 0}]]

    compute_insert_size_threshold()

    for read in reads:
        if not read.is_unmapped and read.positions:
            start, end = true_position(read)
            chromosome = read.reference_name
            # place flag if read and its mate have the same orientation.
            flags, isbuildingflags, all_flags = flag_sameorientation(read, flags, isbuildingflags, all_flags,
                                                                     chromosome, start)

            # place flag if read has high insert size.
            flags, isbuildingflags, all_flags = flag_high_isize(read, flags, isbuildingflags, all_flags, chromosome,
                                                                start)
            # place flag if read has unmapped mate.
            flags, isbuildingflags, all_flags = flag_unmapped_mate(read, flags, isbuildingflags, all_flags, chromosome,
                                                                   start)

            # place flag if read has ultra high insert size.
            flags, isbuildingflags, all_flags = flag_ultra_high_isize(read, flags, isbuildingflags, all_flags,
                                                                      chromosome,  start)

            # place flag if the mate of the read is mapped to a different chromosome.
            flags, isbuildingflags, all_flags = flag_interchromosomal_mate(read, flags, isbuildingflags, all_flags,
                                                                     chromosome, start)

            # place flag if the read and its mate face away from each other
            flags, isbuildingflags, all_flags = flag_facaway(read, flags, isbuildingflags, all_flags,
                                                                           chromosome, start)

            flags = update_total(flags, isbuildingflags)

        elif read.is_unmapped:
            read_data[1] += 1
        elif len(read.positions) == 0:
            read_data[2] += 1
        read_data[0] += 1

    return all_flags, read_data


def compute_insert_size_threshold():
    """ The compute_insert_size_threshold computes the threshold of the high_insert_size flag if the user specified it.
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')
    reads = bamfile.fetch()

    insert_sizes = []
    max_reads = 10000

    global high_threshold
    global ultra_high_threshold

    for read in reads:
        if read.is_proper_pair:
            insert_sizes.append(abs(read.isize))

        if len(insert_sizes) == max_reads:
            break

    insert_sizes = np.array(insert_sizes)
    if args.high_insert_size == -1:
        high_threshold = int(np.percentile(insert_sizes, 99.5))
    else:
        high_threshold = args.high_insert_size


def flag_sameorientation(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_sameorientation function checks if the current read should be added to a same_orientation flag or
    start creating a same_orientation flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    if issameorientation(read):
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 0)

    elif isbuildingflags[0] and start > flags[0][2]:
        percentage = round(flags[0][3]['count'] / flags[0][3]['total'], 2)
        if flags[0][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[0])
        flags[0] = [chromosome, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}]
        isbuildingflags[0] = False

    return flags, isbuildingflags, all_flags


def flag_high_isize(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_high_isize function checks if the current read should be added to a high_insert_size flag or start
    creating a high_isize_flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    insert_size = abs(read.isize)

    if high_threshold < insert_size:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 1)
        flags[1][3]['lengths'].append(insert_size)

    elif isbuildingflags[1] and start > flags[1][2]:
        percentage = round(flags[1][3]['count'] / flags[1][3]['total'], 2)
        if flags[1][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[1])
        flags[1] = [chromosome, None, None, {'type': 'high_insert_size', 'count': 0, 'total': 0, 'lengths': []}]
        isbuildingflags[1] = False

    return flags, isbuildingflags, all_flags


def flag_ultra_high_isize(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_ultra_high_isize function checks if the current read should be added to a Ultra_high_insert_size flag
    or start creating an ultra_high_isize_flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    insert_size = abs(read.isize)

    if insert_size > args.ultra_high_insert_size:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 3)
        flags[3][3]['lengths'].append(insert_size)

    elif isbuildingflags[3] and start > flags[3][2]:
        percentage = round(flags[3][3]['count'] / flags[3][3]['total'], 2)
        if flags[3][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[3])
        flags[3] = [chromosome, None, None, {'type': 'ultra_high_insert_size', 'count': 0, 'total': 0, 'lengths': []}]
        isbuildingflags[3] = False

    return flags, isbuildingflags, all_flags


def flag_unmapped_mate(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_unmapped_mate function checks if the current read should be added to the unmapped_mate flag or start
    creating a unmapped_mate flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    if read.mate_is_unmapped:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 2)

    elif isbuildingflags[2] and start > flags[2][2]:
        percentage = round(flags[2][3]['count'] / flags[2][3]['total'], 2)
        if flags[2][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[2])
        flags[2] = [chromosome, None, None, {'type': 'unmapped_mate', 'count': 0, 'total': 0}]
        isbuildingflags[2] = False

    return flags, isbuildingflags, all_flags


def flag_interchromosomal_mate(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_interchromosomal_mate function checks if the current read should be added to a inter_chromosomal_mate
     flag or start creating a inter_chromosomal flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    if read.reference_name != read.next_reference_name:
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 4)

    elif isbuildingflags[4] and start > flags[4][2]:
        percentage = round(flags[4][3]['count'] / flags[4][3]['total'], 2)
        if flags[4][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[4])
        flags[4] = [chromosome, None, None, {'type': 'inter_chromosomal_pair', 'count': 0, 'total': 0}]
        isbuildingflags[4] = False

    return flags, isbuildingflags, all_flags


def flag_facaway(read, flags, isbuildingflags, all_flags, chromosome, start):
    """ The flag_facaway function checks if the current read should be added to a face_away flag or start creating a
    face_away flag.

    :param read: pysam object containing data of a read.
    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    :param all_flags: A 2d list containing the coordinates of the flags and additional information.
    :param chromosome: The chromosome where the read is mapped.
    :param start: Integer indicating the starting position of the read.

    :return flags: a 2d list containing all the flag information.
    :return isbuildingflags: a list indicating which flags are currently being built.
    :return all_flags: A 2d list containing the coordinates of the flags and additional information.
    """
    if is_facaway(read):
        flags, isbuildingflags = generate_flag(read, flags, isbuildingflags, 5)

    elif isbuildingflags[5] and start > flags[5][2]:
        percentage = round(flags[5][3]['count'] / flags[5][3]['total'], 2)
        if flags[5][3]['count'] > args.threshold and percentage > args.minpercentage:
            all_flags.append(flags[5])
        flags[5] = [chromosome, None, None, {'type': 'face_away', 'count': 0, 'total': 0}]
        isbuildingflags[5] = False

    return flags, isbuildingflags, all_flags


def update_total(flags, isbuildingflags):
    """ The update_total function iterates over all the flag types and increments the total number of reads it has
    encountered by 1.

    :param flags: a 2d list containing all the flag information.
    :param isbuildingflags: a list indicating which flags are currently being built.
    return flags: a 2d list containing all the flag information.
    """
    for index in range(0, len(flags)):
        if isbuildingflags[index]:
            flags[index][3]['total'] += 1

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
        flags[flagindex][0] = read.reference_name
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
    start = read.positions[0] - calculate_overshoot(cigar)
    end = read.positions[-1] + calculate_overshoot(cigar[::-1])

    return start, end


def calculate_overshoot(cigar):
    """ The calculate overshoot function calculates the number of basepairs that have not been mapped but are part of
    the read.

    :param cigar: a list containing tuples representing the cigar string.
    :return overshoot: an integer indicating the number of basepairs in the read before it is mapped.
    """
    overshoot = 0
    for element in cigar:
        if element[0]:
            overshoot += element[1]
        else:
            break

    return overshoot


def issameorientation(read):
    """ The issameorientation function returns a bool returning true if a pair of reads have the same orientation and
    the mate is on the same chromosome. It will return false if this is not the case

    :param read: Pysam object containing data of a read.
    :return bool: A boolean returning True if the reads of a pair have the same orientation.
    """
    if (read.is_paired and not read.mate_is_unmapped and read.is_reverse == read.mate_is_reverse and
        read.reference_name == read.next_reference_name):
        return True

    return False


def is_facaway(read):
    """ The is_facaway function receives a read and returns true if the read and its paired mate face away from each
    other.

    :param read: Pysam object containing data of a read.
    :return bool: A boolean returning True if the read and its mate face away from each other.
    """
    if not issameorientation(read):
        if not read.is_reverse:
            if read.template_length < 0:
                return True
        else:
            if read.template_length > 0:
                return True

    return False


def write_bedfile(flags):
    """ The write_bedfile function writes a file in BED format that can be loaded in igv and visualises the read data.

    :param flags: a 2d list containing all the flag information.
    """
    with open(args.output + f'/{args.name}.bed', 'w') as bedfile:
        bedfile.write('track name=Flags description="Flags regions of interest." db=hg19 gffTags=on itemRGB="On"\n')

    for flag in flags:
        percentage = round(flag[3]['count'] / flag[3]['total'], 2)

        region = f"{flag[0]}\t{flag[1]}\t{flag[2]}"
        description = f"Name={flag[3]['type']};Readcount={flag[3]['count']};Total_reads={flag[3]['total']};" \
                      f"Percentage={percentage}"

        if flag[3]['type'] == 'high_insert_size' or flag[3]['type'] == 'ultra_high_insert_size':
            lengths = flag[3]['lengths']
            description += f";Avg_insert_size={round(mean(lengths))};Med_insert_size={median(lengths)};" \
                           f"Lower_limit={min(lengths)};Upper_limit={max(lengths)}"

        if flag[3]['type'] == 'same_orientation':
            rgb = f"0,192,199"

        elif flag[3]['type'] == 'high_insert_size':
            description += f";Threshold={high_threshold}bp"
            rgb = f"232,135,26"

        elif flag[3]['type'] == 'unmapped_mate':
            rgb = f"218,52,144"

        elif flag[3]['type'] == 'ultra_high_insert_size':
            rgb = '71,226,111'

        elif flag[3]['type'] == 'inter_chromosomal_pair':
            rgb = '87,61,219'

        elif flag[3]['type'] == 'face_away':
            rgb = '147,133,255'

        with open(args.output + f'/{args.name}.bed', 'a') as bedfile:
            bedfile.write(f"{region}\t{description}\t0\t.\t{flag[1]}\t{flag[2]}\t{rgb}\n")


def sort_flags(flags):
    """ The sort_flags function sorts the flags on starting position using insertionsort.

    :param flags: a 2d list containing all the flag information.
    :return flags: a 2d list containing all the flag information.
    """
    for i in range(1, len(flags)):
        key = flags[i][1]

        j = i - 1
        while j >= 0 and key < flags[j][1] and flags[i][0] == flags[j][0]:
            temp = flags[j+1]
            flags[j+1] = flags[j]
            flags[j] = temp
            j -= 1
        flags[j + 1][1] = key

    return flags


def write_logfile(read_data):
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}/Flag_placer.py\nScript finished at: {current_time} {current_day}\n' \
           f'{"-"*40}Read data{"-"*40}\nTotal reads: {read_data[0]}\nUnmapped reads: {read_data[1]}\n' \
           f'Reads without matches: {read_data[2]}\n{"-"*40}Parameters{"-"*40}\nRegion: {args.region}\n' \
           f'Bamfile: {args.bam}\nOutput_folder: {args.output}\nRead threshold: {args.threshold}\n' \
           f'Insert size threshold: {high_threshold}\nUltra high insert size threshold: {args.ultra_high_insert_size}' \
           f'\nMinimal percentage: {args.minpercentage}'

    with open(args.output + f'/{args.name}_log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    """ The main calls the other functions in the right order."""

    reads = fetch_reads()  # get all reads at the regions from the vcf call.

    flags, read_data = place_flags(reads)  # retrieve data from reads.

    sorted_flags = sort_flags(flags)  # sort flags on start position.

    write_bedfile(sorted_flags)  # write the result in a BED file.

    if args.log:
        write_logfile(read_data)  # write logfile with parameters
