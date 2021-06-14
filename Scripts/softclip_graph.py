from datetime import datetime, date
import pysam
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--region', '-r', required=True, type=str,
                    help='String specifying the region in format: "chr#:start-stop". use chr# for whole chromosome.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=False, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--name', '-n', required=False, default='output', type=str,
                    help='Name for the project. This is the name of the output file.')
args = parser.parse_args()


def fetch_reads():
    """ The fetch_reads function fetches the reads from the bam file

    :return reads: Pysam object containing read information.
    :return pileup: Pysam object containing read information of all positions in the region.
    """
    global bamfile
    bamfile = pysam.AlignmentFile(args.bam, 'rb')

    if ':' in args.region and '-' in args.region:
        chromosome, start, end = re.split(':|-', args.region)
        chromosome = chromosome.replace('chr', '')

        reads = bamfile.fetch(chromosome, int(start), int(end))
        pileup = bamfile.pileup(chromosome, int(start), int(end))

    else:
        chromosome = args.region.replace('chr', '')
        reads = bamfile.fetch(chromosome)
        pileup = bamfile.pileup(chromosome)

    return reads, pileup


def get_softclipdata(reads):
    """ The get_softclipdata function is a function that retrieves the position of all softclip bases. It also keeps
    track of the number of the total number of reads and the number of reads it could not interpret.

    :param reads: A pysam object containing all the reads with their information.
    :return softclipdata: A dictionary with positions and the number of "normal" bases and softclip bases.
    :return read_data: a list containing the total number of reads and the reads it could not interpret.
    """
    softclipdata = {}
    read_data = [0, 0, 0]

    for read in reads:
        if not read.is_unmapped and len(read.positions) != 0:
            if 'S' in read.cigarstring:
                softclip_pos = get_softclip_positions(read)
                softclipdata = update_softclipdata(softclipdata, softclip_pos)

        elif read.is_unmapped:
            read_data[1] += 1
        elif len(read.positions) == 0:
            read_data[2] += 1
        read_data[0] += 1

    return softclipdata, read_data


def get_softclip_positions(read):
    """ The get_softclip_positions function receives a read and returns a list with the positions of softclipped bases.

    :param read: Pysam object containing read information.
    :return softclip_pos: A list containing the positions of softclipped bases.
    """
    softclip_pos = []

    cursor = true_start(read.cigar, read.positions[0])
    for element in read.cigartuples:
        if element[0] == 4:
            softclip_pos += [x for x in range(cursor, cursor + element[1])]

        cursor += element[1]

    return softclip_pos


def true_start(cigar, matchstart):
    """ The true_start function receives the cigar string and the starting position of the first match in a read. It
    returns the start position of the read including unmapped parts.

    :param cigar: a list containing tuples representing the cigar string.
    :param matchstart: an integer representing the start of the first mapped base of a read.
    :return read_start: an integer representing the starting position of the read including softclips.
    """
    overshoot = 0

    for element in cigar:
        if element[0] != 0:
            overshoot += element[1]
        else:
            break

    read_start = matchstart - overshoot

    return read_start


def update_softclipdata(softclipdata, softclip_pos):
    """ The update_softclipdata function updates the softclipdata dictionary when a new softclipped base is found.

    :param softclipdata: A dictionary with positions and the number of softclip bases.
    :param softclip_pos: An integer representing the position of the softclip base on the reference genome.
    :return softclipdata: A dictionary with positions and the number of softclip bases.
    """
    for pos in softclip_pos:
        if pos in softclipdata:
            softclipdata[pos][0] += 1
        else:
            softclipdata.update({pos: [1, 0]})

    return softclipdata


def get_coverage(ref, pos):
    """ The get_coverage function retrieves the number of aligned reads of a position.

    :param ref: a string representing the chromosome.
    :param pos: an integer representing the position of the reference genome.
    :return coverage: an integer representing the coverage of the given position.
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')
    coverage = 0
    for pileupcolumn in bamfile.pileup(ref, pos, pos + 1):
        coverage = pileupcolumn.nsegments
        break

    return coverage


def add_coverage(softclipdata, pileup):
    """ The add coverage function adds the coverage to the softclipdata dictionary.

    :param softclipdata: A dictionary with positions and the number of softclip bases.
    :param pileup: Pysam object containing read information of all positions in the region.
    :return softclipdata: A dictionary with positions with the number of softclip bases and coverage.
    """
    for pileupcolumn in pileup:
        if pileupcolumn.pos in softclipdata:
            softclipdata[pileupcolumn.pos][1] = pileupcolumn.nsegments

    return softclipdata


def to_graphdata(softclipdata):
    """ The to_graphdata function receives the softclipdata dictionary and changes it to a 2d list datastructure that is
    sorted on starting position and compressed by merging positions with the same ratio.

    :param softclipdata: A dictionary with positions with the number of softclip bases and coverage.
    :return compressed_graphdata: A 2d list containing the coordinates and the percentage of sofclipped bases.
    """
    graphdata = []
    chromosome = args.region.split(':')[0]

    for pos in softclipdata:
        percentage = round(softclipdata[pos][0] / (softclipdata[pos][0] + softclipdata[pos][1]), 2)
        graphdata.append([chromosome, pos, pos, percentage])

    graphdata = sorted(graphdata, key=lambda x: x[1])

    compressed_graphdata = [graphdata[0]]

    for row in graphdata[1:]:
        if row[3] == compressed_graphdata[-1][3] and row[2] == compressed_graphdata[-1][2] + 1:
            compressed_graphdata[-1][2] = row[2]
        else:
            compressed_graphdata.append(row)

    return compressed_graphdata


def write_bedgraph_file(graphdata):
    """ The write_bedgraph_file function receives the heatmapdata and writes a BedGraph file.

    :param heatmapdata: A 2d list containing the coordinates and the percentage of sofclipped bases.
    """
    with open(args.output + f'/{args.name}.BedGraph', 'w') as bedfile:
        bedfile.write('track type=bedGraph name=Softclip_graph description="Softclip graph" color=220,20,60 '
                      'graphType=bar alwaysZero=off autoScale=on\n')

    with open(args.output + f'/{args.name}.BedGraph', 'a') as bedfile:
        for datapoint in graphdata:
            bedfile.write(f"{datapoint[0]}\t{datapoint[1]}\t{datapoint[2]}\t{datapoint[3]}\n")


def write_logfile(read_data):
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}/softclip_graph.py\nScript finished at: {current_time} {current_day}\n' \
           f'{"-" * 40}Read data{"-" * 40}\nTotal reads: {read_data[0]}\nUnmapped reads: {read_data[1]}\n' \
           f'Reads without matches: {read_data[2]}\n{"-" * 40}Parameters{"-" * 40}\nRegion: {args.region}\n' \
           f'Bamfile: {args.bam}\nOutput_folder: {args.output}\n'

    with open(args.output + f'/{args.name}_BedGraph_log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    reads, pileup = fetch_reads()  # fetch reads and pileup engine

    softclipdata, read_data = get_softclipdata(reads)  # retreive number of softclips and read information

    softclipdata = add_coverage(softclipdata, pileup)  # add coverage to softclipdata to calculate ratio

    compressed_graphdata = to_graphdata(softclipdata)  # change datastructure of softclipdata

    write_bedgraph_file(compressed_graphdata)  # write to BedGraph formatted file

    if args.log:
        write_logfile(read_data)  # write logfile
