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

args = parser.parse_args()


def fetch_reads():
    """ The fetch_reads function fetches the reads from the bam file

    :return reads: Pysam object containing read information
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')

    reads = bamfile.fetch('3')

    return reads


def place_flags(reads):
    """ The place_flags function gets the coordinates of interesting regions and returns the coordinates of the flags
    with their information.

    :param reads: Pysam object containing read information
    :return flags: A 2d list containing the coordinates of the flags and additional information.
    """
    flags = []

    isbuildingflag = False
    flag = ['3', None, None, 1]

    for read in reads:
        if not read.is_unmapped:
            if issameorientation(read) and not isbuildingflag:
                flag[1] = read.positions[0]
                flag[2] = read.positions[-1]
                isbuildingflag = True

            elif issameorientation(read) and isbuildingflag:
                flag[2] = read.positions[-1]
                flag[3] += 1

            elif not issameorientation(read) and isbuildingflag:
                if read.positions[0] > flag[2]:
                    if flag[3] > args.threshold:
                        flags.append(flag)
                    flag = ['3', None, None, 1]
                    isbuildingflag = False

    return flags


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
    text = 'track name=same_direction_reads description="Region_Summary." db=hg19 gffTags=on\n'

    for flag in flags:
        region = f"{flag[0]}\t{flag[1]}\t{flag[2]}"
        text += f"{region}\tName=Same_facing_reads;Readcount={flag[3]}\n"

    with open(args.output + '/outputheatmap.bed', 'w') as bedfile:
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
