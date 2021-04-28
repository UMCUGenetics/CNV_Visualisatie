import pysam
from pysam import VariantFile
from datetime import datetime, date
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bamfile', '-bf', required=True, type=str, help='Path to bam file.')
parser.add_argument('--vcf_file', '-vcf', required=True, type=str, help='Path to vcf file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--logfile', '-l', required=False, default=True, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--capture_region', '-cr', required=False, default=3000, type=int,
                    help='Number of basepairs that should be included around CNV call.')
parser.add_argument('--high_insert_size', '-hi', required=False, default=500, type=int,
                    help='Length of insert size to be classified as high.')

args = parser.parse_args()


def vcf_calls():
    """ The vcf_calls function opens a vcf file and retrieves al the reqions called in the vcf file.

    :return regions: A list of coordinates specified in the vcf file.
    """
    vcf_file = VariantFile(args.vcf_file)
    regions = []

    for call in vcf_file:
        splitcall = str(call).split('\t')
        chr = splitcall[0]
        start = str(int(splitcall[1]) - args.capture_region)
        stop = str(int(splitcall[7].split(';')[2].split('=')[1]) + args.capture_region)
        region = [chr, start, stop]
        regions.append(region)

    return regions


def fetch_reads(regions):
    """ The fetch_reads function fetches the reads and matches them with their pairs from the regions specified in the
    regions list.

    :param regions: A list containing the chromosome and coordinates for every region.
    :return all_reads: A list containing the number of reads in the specified region
    """
    samfile = pysam.AlignmentFile(args.bamfile, 'rb')
    all_reads = []

    for loc in regions:
        reads = samfile.fetch(str(loc[0]), int(loc[1]), int(loc[2]))

        regional_reads = []
        for read in reads:
            regional_reads.append(read)

        all_reads.append(regional_reads)

    samfile.close()

    return all_reads


def get_read_data(all_reads):
    """ The get_read_data function collects the information from the reads and returns it.

    :param all_reads: a 2d list containing all the reads for every region.
    :return data: a 2d list containing all the data from the reads for each region.
    """
    data = []

    for region in all_reads:
        paired_reads, unmapped_mate, duplicate_pairs, high_isize, facaway, same_orientation, proper_pair = get_read_stats(region)

        regiondata = [paired_reads, unmapped_mate, duplicate_pairs, high_isize, facaway, same_orientation, proper_pair]
        data.append(regiondata)

    return data


def get_read_stats(region):
    """ The get_count_stats function receives all reads in every region and counts the total various features of the
    reads.

    :param region: a list containing the reads of a region.
    :return readcount: an integer representing the total number of reads in the region.
    """
    already_done = {}

    paired_reads = 0
    proper_pair = 0
    unmapped_mate = 0
    duplicate_pairs = 0
    high_isize = 0
    facaway = 0
    same_orientation = 0

    for read in region:
        if read.is_paired and read.qname not in already_done:
            paired_reads += 1

        if read.is_proper_pair and read.qname not in already_done:
            proper_pair += 1

        if read.mate_is_unmapped:
            unmapped_mate += 1

        if read.is_duplicate and read.qname not in already_done:
            duplicate_pairs += 1

        if isvalidread(read, already_done):
            if not -args.high_insert_size <= read.isize <= args.high_insert_size:
                high_isize += 1

        if isfacaway(read) and read.qname not in already_done:
            facaway += 1

        if issameorientation(read) and read.qname not in already_done:
            same_orientation += 1

        already_done.update({read.qname: None})

    return paired_reads, unmapped_mate, duplicate_pairs, high_isize, facaway, same_orientation, proper_pair


def isfacaway(read):
    """ The isfacaway function returns True if the reads in a pair are faced away from each other and false if they are
    not.

    :param read: Pysam object containing data of a read.
    :return bool: A bool telling if a read pair face away from each other.
    """
    if read.is_paired and not read.mate_is_unmapped:
        if read.is_read1:
            if read.is_reverse and not read.mate_is_reverse:
                return True
            else:
                return False

        else:
            if not read.is_reverse and read.mate_is_reverse:
                return True
            else:
                return False

    else:
        return False


def isvalidread(read, already_done):
    """ Function that returns True if read is valid to be used for the high insert size calculation.

    :param read: Pysam object containing data of a read.
    :param already_done: Dictionary containing read names.
    :return boool: A boolean returning true if read is valid.
    """
    if (read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and not read.is_duplicate and
            read.qname not in already_done):
        return True
    else:
        return False


def issameorientation(read):
    """ The issameorientation function returns a bool returning true if a pair of reads have the same orientation and
    false if they have an oposite orientation.

    :param read: Pysam object containing data of a read.
    :return bool: A boolean returning True if the reads of a pair have the same orientation.
    """
    if read.is_paired and not read.mate_is_unmapped:
        if read.is_reverse and read.mate_is_reverse:
            return True

        elif not read.is_reverse and not read.mate_is_reverse:
            return True

        else:
            return False


def write_bedfile(regions, read_data):
    """ The write_bedfile function writes a file in BED format that can be loaded in igv and visualises the read data.

    :param regions: a list of coordinates specified in the vcf file.
    :param read_data: a 2d list containing the read data of every region.
    """
    text = 'track name=CNV_information description="Region_Summary." db=hg19 gffTags=on\n'

    for index in range(0, len(regions)):
        region = f"{regions[index][0]}\t{regions[index][1]}\t{regions[index][2]}"
        paired_reads = read_data[index][0]
        unmapped_mate = read_data[index][1]
        duplicate_pairs = read_data[index][2]
        high_insize = read_data[index][3]
        facaway = read_data[index][4]
        same_orientation = read_data[index][5]
        proper_pair = read_data[index][6]

        text += f"{region}\tName=Read_information;Paired_reads={paired_reads};Proper_pairs={proper_pair};Unmapped_mate=" \
                f"{unmapped_mate};Duplicate_pair={duplicate_pairs};High_insert_size={high_insize};Facaway={facaway};" \
                f"Same_orientation={same_orientation}\n"

    with open(args.output + '/output.bed', 'w') as igvfile:
        igvfile.write(text)


def write_logfile():
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}\nScript finished at: {current_time} {current_day}\n{"-"*80}\n' \
           f'Parameters:\nBamfile: {args.bamfile}\nVCF file: {args.vcf_file}\nOutput_folder: {args.output}\n' \
           f'Capture region: {args.capture_region}\nHigh insert size threshold: {args.high_insert_size}'

    with open(args.output + '/log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    """ The main calls the other functions in the right order."""

    regions = vcf_calls()  # get coordinates of the regions specified in vcf file.

    all_reads = fetch_reads(regions)  # get all reads at the regions from the vcf call.

    read_data = get_read_data(all_reads)  # retrieve data from reads.

    write_bedfile(regions, read_data)  # write the result in a BED file.

    if args.logfile:
        write_logfile()  # write logfile with parameters
