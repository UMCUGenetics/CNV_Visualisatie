import os
import argparse
import time
from datetime import datetime, date
from multiprocessing import Pool


parser = argparse.ArgumentParser()
parser.add_argument('--bam', '-b', required=True, type=str, help='Path to bam file.')
parser.add_argument('--output', '-o', required=False, type=str, help='Path to output folder.')
parser.add_argument('--name', '-n', required=False, default='output', type=str,
                    help='Name for the project. This is the name of the output file.')
parser.add_argument('--cores', '-c', required=False, default=1, type=int, help='Number of cores that should be used.')

args = parser.parse_args()


def read_settings():
    """ The read_settings function reads the settings file and creates a global dictionary containing all the settings.
    """
    global settings
    settings = {}

    with open('../Settings.txt', 'r') as file:
        text = file.readlines()

    for row in text:
        if not row.startswith('#'):
            if row != '\n':
                row.replace('\n', '')
                splitrow = row.replace('\n', '').split('=')
                settings.update({splitrow[0]: splitrow[1]})


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
        os.system(f'python3 Softclip_graph.py -b "{args.bam}"'
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
    if extension == 'bed':
        type = 'flags'
    else:
        type = 'softclip'

    try:
        with open(f"{args.output}/{args.name}_{type}_{chromosome}.{extension}", 'r') as bedfile:
            text = bedfile.readlines()[1:]
        for line in text:
            bedfile_text += line

        os.remove(f"{args.output}/{args.name}_{type}_{chromosome}.{extension}")

    except FileNotFoundError:
        print(f"WARNING could not merge file: '{args.output}/{args.name}_{type}_{chromosome}.{extension}'")

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


def get_log_data(chromosomes):
    total_reads = 0
    total_unmapped = 0
    total_reads_without_matches = 0

    for chromosome in chromosomes:
        with open(f"{args.output}/{args.name}_flags_{chromosome}_log.txt", 'r') as bedfile:
            text = bedfile.readlines()

        total_reads += int(text[3].split(': ')[1])
        total_unmapped += int(text[4].split(': ')[1])
        total_reads_without_matches += int(text[5].split(': ')[1])

        os.remove(f"{args.output}/{args.name}_{chromosome}_bed_log.txt")
        os.remove(f"{args.output}/{args.name}_{chromosome}_BedGraph_log.txt")

        return total_reads, total_unmapped, total_reads_without_matches


def merge_logfiles(chromosomes):
    total_reads, total_unmapped, total_reads_without_matches = get_log_data(chromosomes)

    current_path = os.getcwd()

    current_time = datetime.now().strftime("%H:%M:%S")
    current_day = date.today().strftime("%d/%m/%Y")

    text = f'Logfile created by: {current_path}/Start_job.py\nScript finished at: {current_time} {current_day}\n' \
           f'{"-" * 40}Read data{"-" * 40}\nTotal reads: {total_reads}\nUnmapped reads: {total_unmapped}\n' \
           f'Reads without matches: {total_reads_without_matches}\n{"-" * 40}Parameters{"-" * 40}\nRegion: {args.region}' \
           f'\nBamfile: {args.bam}\nOutput_folder: {args.output}\nCores: {args.cores}\n{"-" * 40}Settings{"-" * 40}' \
           f'\nhigh_insert_size={settings["high_insert_size"]}\nultra_high_insert_size=' \
           f'{settings["ultra_high_insert_size"]}\nMinPercentage_same_orientation=' \
           f'{settings["MinPercentage_same_orientation"]}\nMinPercentage_high_insert_size=' \
           f'{settings["MinPercentage_high_insert_size"]}\nMinPercentage_unmapped_mate=' \
           f'{settings["MinPercentage_unmapped_mate"]}\nMinPercentage_ultra_high_insert_size=' \
           f'{settings["MinPercentage_ultra_high_insert_size"]}\nMinPercentage_inter_chromosomal_pair=' \
           f'{settings["MinPercentage_inter_chromosomal_pair"]}\nMinPercentage_face_away=' \
           f'{settings["MinPercentage_face_away"]}\nMinReadCount_same_orientation=' \
           f'{settings["MinReadCount_same_orientation"]}\nMinReadCount_high_insert_size=' \
           f'{settings["MinReadCount_high_insert_size"]}\nMinReadCount_unmapped_mate=' \
           f'{settings["MinReadCount_unmapped_mate"]}\nMinReadCount_ultra_high_insert_size=' \
           f'{settings["MinReadCount_ultra_high_insert_size"]}\nMinReadCount_inter_chromosomal_pair=' \
           f'{settings["MinReadCount_inter_chromosomal_pair"]}\nMinReadCount_face_away=' \
           f'{settings["MinReadCount_face_away"]}\nMinCoverage={settings["MinCoverage"]}\n'

    with open(args.output + f'/{args.name}_log.txt', 'w') as logfile:
        logfile.write(text)


if __name__ == '__main__':
    start = time.time()  # Keep track of time.
    chromosomes = list(range(1, 23)) + ['X', 'Y']  # create a list of all chromosomes.
    read_settings()

    bedfile_handle(chromosomes, 'flags')  # create bed files for each chromosome.

    merge_bedfiles(chromosomes, 'bed')  # merge bed files created by bedfile_handle.

    bedfile_handle(chromosomes, 'softclips')  # create BedGraph files for each chromosome.

    merge_bedfiles(chromosomes, 'BedGraph')  # merge BedGraph files created by bedgraph_handle.

    if settings['log'] == 'True':
        merge_logfiles(chromosomes)

    end = time.time()  # stop tracking time

    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"elapsed time: {round(hours)}:{round(minutes)}:{round(seconds)}")  # print elapsed time.
