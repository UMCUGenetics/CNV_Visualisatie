import os
import argparse
import time
from multiprocessing import Pool


parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--log', '-l', required=False, default=False, type=bool,
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
parser.add_argument('--name', '-n', required=False, default='output', type=str, help='Name for the project. This is'
                                                                                         'the name of the output file.')
parser.add_argument('--cores', '-c', required=False, default=1, type=int, help='Number of cores that should be used.')

args = parser.parse_args()


def write_bedfile(chromosome):
    """ The write_bedfile function runs the Flag_placer.py script for the given chromosome with the given arguments.

    :param chromosome: Int or Str specifying the chromosome.
    """
    os.system(f'python3 Flag_placer.py -b "{args.bam}"'
              f' -o "{args.output}"'
              f' -r "chr{chromosome}:0-0"'
              f' -t {args.threshold}'
              f' -mp {args.minpercentage}'
              f' -n "{args.name}_{chromosome}"')


def write_bedgraphfile(chromosome):
    """ The write_bedgraphfile function runs the softclipp_graph.py script for the given chromosome with the given
    arguments

    :param chromosome: Int or Str specifying the chromosome.
    """
    os.system(f'python3 softclip_graph.py -b "{args.bam}"'
              f' -o "{args.output}"'
              f' -r "chr{chromosome}:0-0"'
              f' -n "{args.name}_{chromosome}"')


def merge_bedfiles(chromosomes, extension):
    """ The merge_bedfiles function combines all the bedfiles or BedGraph files into one large file.

    :param chromosomes: Int or Str specifying the chromosome.
    :param extension: Str representing the file extension.
    """
    bedfile_text = 'track name=Flags description="Flags regions of interest." db=hg19 gffTags=on itemRGB="On"\n'

    for chromosome in chromosomes:
        with open(f"{args.output}/{args.name}_{chromosome}.{extension}", 'r') as bedfile:
            text = bedfile.readlines()[1:]
        for line in text:
            bedfile_text += line

        os.remove(f"{args.output}/{args.name}_{chromosome}.{extension}")

    with open(f"{args.output}/{args.name}.{extension}", 'w') as output:
        output.write(bedfile_text)


def bedfile_handle(chromosomes):
    with Pool(args.cores) as p:
        p.map(write_bedfile, chromosomes)


def bedgraph_handle(chromosomes):
    with Pool(args.cores) as p:
        p.map(write_bedgraphfile, chromosomes)


if __name__ == '__main__':
    start = time.time()
    chromosomes = list(range(1, 23)) + ['X', 'Y']

    bedfile_handle(chromosomes)

    merge_bedfiles(chromosomes, 'bed')

    bedgraph_handle(chromosomes)

    merge_bedfiles(chromosomes, 'BedGraph')

    end = time.time()

    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"elapsed time: {hours}:{minutes}:{seconds}")
