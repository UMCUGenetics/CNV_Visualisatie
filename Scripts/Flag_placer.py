import pysam
from datetime import datetime, date
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=True, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--threshold', '-t', required=False, default=0, type=int, help='int specifying the minimun amount'
                                                                                    'of reads before a flag is made')
parser.add_argument('--region', '-r', required=False, default='', type=str, help='String specifying the region in '
                                                                                 'format: "chr#:start-stop". use '
                                                                                 'chr#:0-0 for whole chromosome.')

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
    :return flags: A 2d list containing the coordinates of the flags and additional information.
    """
    all_flags = []
    chromosome = args.region.split(':')[0].replace('chr', '')

    isbuildingflags = [False, False]
    flags = [[chromosome, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}],
             [chromosome, None, None, {'type': str, 'count': 0, 'total': 0}]]

    for read in reads:
        if not read.is_unmapped:
            if issameorientation(read):
                flags, isbuildingflags = same_orientation_flag(read, flags, isbuildingflags)

            elif not issameorientation(read) and isbuildingflags[0] and read.positions[0] > flags[0][2]:
                if flags[0][3]['count'] > args.threshold:
                    all_flags.append(flags[0])
                flags[0] = [chromosome, None, None, {'type': 'same_orientation', 'count': 0, 'total': 0}]
                isbuildingflags[0] = False

            flags = update_total(flags)

    return all_flags


def update_total(flags):

    for flag in flags:
        total = flag[3]['total']
        total += 1
        flag[3].update({'total': total})

    return flags


def same_orientation_flag(read, flags, isbuildingflags):
    if not isbuildingflags[0]:
        flags[0][1] = read.positions[0]
        flags[0][2] = read.positions[-1]
        isbuildingflags[0] = True

    if isbuildingflags[0]:
        flags[0][2] = read.positions[-1]

    counter = flags[0][3]['count']
    counter += 1
    flags[0][3].update({'count': counter})

    return flags, isbuildingflags


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


def write_bedfile(flags):
    """ The write_bedfile function writes a file in BED format that can be loaded in igv and visualises the read data.

    :param regions: a list of coordinates specified in the vcf file.
    :param read_data: a 2d list containing the read data of every region.
    """
    # TODO: add color scaling to rgb
    text = 'track name=same_direction_reads description="Region_Summary." db=hg19 gffTags=on itemRGB="On"\n'

    for flag in flags:
        percentage = round(flag[3]['count'] / flag[3]['total'], 2)

        region = f"{flag[0]}\t{flag[1]}\t{flag[2]}"
        description = f"Name={flag[3]['type']};Readcount={flag[3]['count']};Total_reads={flag[3]['total']};" \
                      f"Percentage={percentage}"
        if flag[3]['type'] == 'same_orientation':
            rgb = f"150,200,150"
        else:
            rgb = '0,0,0'

        text += f"{region}\t{description}\t0\t.\t{flag[1]}\t{flag[2]}\t{rgb}\n"

    with open(args.output + '/output_flags.bed', 'w') as bedfile:
        bedfile.write(text)


def write_logfile():
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}\nScript finished at: {current_time} {current_day}\n{"-"*80}\n' \
           f'Parameters:\nBamfile: {args.bam}\nOutput_folder: {args.output}\n' \
           f'Read threshold: {args.threshold}\n'

    with open(args.output + '/log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    """ The main calls the other functions in the right order."""

    reads = fetch_reads()  # get all reads at the regions from the vcf call.

    flags = place_flags(reads)  # retrieve data from reads.

    write_bedfile(flags)  # write the result in a BED file.

    if args.log:
        write_logfile()  # write logfile with parameters
