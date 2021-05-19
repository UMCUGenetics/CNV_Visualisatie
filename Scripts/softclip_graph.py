from datetime import datetime, date
import pysam
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=False, type=bool,
                    help='Bool specifying if logfile should be made.')
parser.add_argument('--region', '-r', required=False, default='', type=str, help='String specifying the region in '
                                                                                 'format: "chr#:start-stop". use '
                                                                                 'chr#:0-0 for whole chromosome.')
parser.add_argument('--name', '-n', required=False, default='output', type=str, help='Name for the project. This is'
                                                                                     'the name of the output file.')
args = parser.parse_args()


def fetch_reads():
    """ The fetch_reads function fetches the reads from the bam file

    :return reads: Pysam object containing read information
    """
    bamfile = pysam.AlignmentFile(args.bam, 'rb')

    if args.region == 'all':
        reads = bamfile.fetch()

    else:
        chromosome = args.region.split(':')[0].replace('chr', '')
        start = args.region.split(':')[1].split('-')[0]
        end = args.region.split(':')[1].split('-')[1]

        if start == '0' and end == '0':
            reads = bamfile.fetch(chromosome)

        else:
            reads = bamfile.fetch(chromosome, int(start), int(end))

    return reads


def get_softclipdata(reads):
    """ The get_softclipdata function is a function that retrieves the position of all softclip bases. It also keeps
    track of the number of the total number of reads and the number of reads it could not interpret.

    :param reads: A pysam object containing all the reads with their information.
    :return softclipdata: A dictionary with positions and the number of "normal" bases and softclip bases.
    :read_data: a list containing the total number of reads and the reads it could not interpret.
    """
    softclipdata = {}
    read_data = [0, 0, 0]

    for read in reads:
        if not read.is_unmapped and len(read.positions) != 0:
            read_start = true_start(read.cigar, read.positions[0])
            regions = softclip_regions(read_start, read.cigar)

            if len(regions) != 0:
                for region in regions:
                    softclipdata = update_softclipdata(region, softclipdata)

        elif read.is_unmapped:
            read_data[1] += 1
        elif len(read.positions) == 0:
            read_data[2] += 1
        read_data[0] += 1

    return softclipdata, read_data


def update_softclipdata(region, softclipdata):
    """ The update_softclipdata function updates the softclipdata dictionary with new information.

    :param region: A region of a read defined by the cigar string.
    :param softclipdata: A dictionary with positions and the number of "normal" bases and softclip bases.
    :return sofclipdata: A dictionary with positions and the number of "normal" bases and softclip bases.
    """
    if region[2] == 'normal':
        index = 0
    else:
        index = 1

    for pos in range(region[0], region[1]):
        if pos in softclipdata:
            softclipdata[pos][index] += 1
        else:
            posdata = [0, 0]
            posdata[index] += 1
            softclipdata.update({pos: posdata})

    return softclipdata


def softclip_regions(read_start, cigar):
    """ The softclip regions function iterates over the cigar string and returns the positions and if they are "normal"
    or a softclipped region.

    :param read_start: An integer representing the start of the read (including unmapped bases)
    :param cigar: a list containing tuples representing the cigar string.
    :return regions: a 2d list containing the regions of a read and specifying if they are "normal" or softclipped.
    """
    regions = []
    cursor = read_start
    for element in cigar:
        if element[0] == 4:
            regions.append([cursor, cursor+element[1], 'softclip'])
        else:
            regions.append([cursor, cursor+element[1], 'normal'])
        cursor += element[1]

    return regions


def true_start(cigar, matchstart):
    """ The true_start function receives the cigar string and the starting position of the first match in a read. It
    returns the start position of the read including unmapped parts.

    :param cigar: a list containing tuples representing the cigar string.
    :param matchstart: an integer representing the start of the first mapped base of a read.
    """
    overshoot = 0

    for element in cigar:
        if element[0] != 0:
            overshoot += element[1]
        else:
            break

    read_start = matchstart - overshoot

    return read_start


def heatmapdata(sofclipdata):
    """ The heatmapdata function receives the softclipdata and translates it to a more suitable format for a bedfile.

    :param sofclipdata: A dictionary with positions and the number of "normal" bases and softclip bases.
    :return heatmapdata: A 2d list containing the coordinates and the percentage of sofclipped bases.
    """
    heatmapdata = []
    first = True
    chromosome = args.region.split(':')[0]

    for pos in sofclipdata:
        if sofclipdata[pos][1] != 0:
            ratio = round(sofclipdata[pos][1]/(sofclipdata[pos][0]+sofclipdata[pos][1]), 2)
            if not first:
                if ratio == heatmapdata[-1][3] and heatmapdata[-1][2] == (pos-1):
                    old = heatmapdata[-1]
                    old[2] = pos
                else:
                    heatmapdata.append([chromosome, pos, pos, ratio])
            else:
                heatmapdata.append([chromosome, pos, pos, ratio])
                first = False

    return heatmapdata


def write_bedgraph_file(heatmapdata):
    """ The write_bedgraph_file function receives the heatmapdata and writes a BedGraph file.

    :param heatmapdata: A 2d list containing the coordinates and the percentage of sofclipped bases.
    """
    text = 'track type=bedGraph name=Softclip_graph description="Softclip graph" color=220,20,60 ' \
           'graphType=bar alwaysZero=off\n'
    for datapoint in heatmapdata:
        text += f"{datapoint[0]}\t{datapoint[1]}\t{datapoint[2]}\t{datapoint[3]}\n"

    with open(args.output + f'/{args.name}.BedGraph', 'w') as bedfile:
        bedfile.write(text)


def write_logfile(read_data):
    """ The write logfile function writes a log.txt file in the output folder and writes all the parameters down."""
    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}/softclip_heatmap.py\nScript finished at: {current_time} {current_day}\n' \
           f'{"-"*40}Read data{"-"*40}\nTotal reads: {read_data[0]}\nUnmapped reads: {read_data[1]}\n' \
           f'Reads without matches: {read_data[2]}\n{"-"*40}Parameters{"-"*40}\nRegion: {args.region}\n' \
           f'Bamfile: {args.bam}\nOutput_folder: {args.output}\n'

    with open(args.output + f'/{args.name}_BedGraph_log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    reads = fetch_reads()

    softclipdata, read_data = get_softclipdata(reads)
    heatmapdata = heatmapdata(softclipdata)

    if args.log:
        write_logfile(read_data)

    write_bedgraph_file(heatmapdata)
