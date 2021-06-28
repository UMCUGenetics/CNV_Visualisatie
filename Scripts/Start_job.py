import os
import argparse
import time
from multiprocessing import Pool


parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--name', '-n', required=False, default='output', type=str,
                    help='Name for the project. This is the name of the output file.')
parser.add_argument('--cores', '-c', required=False, default=1, type=int, help='Number of cores that should be used.')

args = parser.parse_args()


def write_bedfile(chromosome):
    """ The write_bedfile function runs the Flag_placer.py script for the given chromosome with the given arguments.

    :param chromosome: Int or Str specifying the chromosome.
    """
    if not os.path.exists(f"{args.output}/{args.name}_{chromosome}.bed"):
        os.system(f'python3 Flag_placer.py -b "{args.bam}"'
                  f' -o "{args.output}"'
                  f' -r "chr{chromosome}"'
                  f' -n "{args.name}_{chromosome}"')


def write_bedgraphfile(chromosome):
    """ The write_bedgraphfile function runs the softclipp_graph.py script for the given chromosome with the given
    arguments

    :param chromosome: Int or Str specifying the chromosome.
    """
    if not os.path.exists(f"{args.output}/{args.name}_{chromosome}.BedGraph"):
        os.system(f'python3 softclip_graph.py -b "{args.bam}"'
                  f' -o "{args.output}"'
                  f' -r "chr{chromosome}"'
                  f' -n "{args.name}_{chromosome}"')


def merge_bedfiles(chromosomes, extension):
    """ The merge_bedfiles function combines all the bedfiles or BedGraph files into one large file.

    :param chromosomes: Int or Str specifying the chromosome.
    :param extension: Str representing the file extension.
    """
    with open(f"{args.output}/{args.name}.{extension}", 'w') as output:
        output.write('track name=Flags description="Flags regions of interest." db=hg19 gffTags=on itemRGB="On"\n')

    for chromosome in chromosomes:
        bedfile_text = get_bed_text(chromosome, extension)

        with open(f"{args.output}/{args.name}.{extension}", 'a') as output:
            output.write(bedfile_text)


def get_bed_text(chromosome, extension):
    """ The update text function receives the already merged bedfile_text and adds the content of the next bedfile to
    the total bedfile_text.

    :param bedfile_text: a string containing the content of the merged bedfile
    :param chromosome: a string resembling a chromosome.
    :param extension: a string resembling the file extension.
    :return bedfile_text: a string containing the text of the bed file.
    """
    bedfile_text = ''
    try:
        with open(f"{args.output}/{args.name}_{chromosome}.{extension}", 'r') as bedfile:
            text = bedfile.readlines()[1:]
        for line in text:
            bedfile_text += line

        os.remove(f"{args.output}/{args.name}_{chromosome}.{extension}")

    except FileNotFoundError:
        print(f"WARNING could not merge file: '{args.output}/{args.name}_{chromosome}.{extension}'")

    return bedfile_text


def bedfile_handle(chromosomes, extension):
    """ The bedfile_handle function divides the chromosomes over the number of cores to multiprocess the Flag_placer.py
    script.

    :param chromosomes: A list of all chromosomes.
    :param extension: A string identifying the script that should be called.
    """
    if extension == 'flags':
        with Pool(args.cores) as p:
            p.map(write_bedfile, chromosomes)

    else:
        with Pool(args.cores) as p:
            p.map(write_bedgraphfile, chromosomes)


if __name__ == '__main__':
    start = time.time()  # Keep track of time.
    chromosomes = list(range(1, 23)) + ['X', 'Y']  # create a list of all chromosomes.

    bedfile_handle(chromosomes, 'flags')  # create bed files for each chromosome.

    merge_bedfiles(chromosomes, 'bed')  # merge bed files created by bedfile_handle.

    bedfile_handle(chromosomes, 'softclips')  # create BedGraph files for each chromosome.

    merge_bedfiles(chromosomes, 'BedGraph')  # merge BedGraph files created by bedgraph_handle.

    end = time.time()  # stop tracking time

    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"elapsed time: {round(hours)}:{round(minutes)}:{round(seconds)}")  # print elapsed time.
