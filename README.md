# CNV_Visualisation

## Description

The CNV_visualisation project is a project about visually detecting and analyzing CNV calls. It contains multiple python
scripts that retrieve additional information from bam files. The CNV_vis.py script creates a BED formatted file that can
be loaded in IGV and gives a summary of the read information around CNV calls made in a vcf file. The flag_placer.py is
a script creates a BED formatted file that highlights potential regions of interest in CNV calling. The softclip_graph.py
script creates a BedGraph formatted file that shows the percentage of softclipped bases for each basepair. The Start_job.py
script can be used to run the flag_placer.py and the softclip_graph.py on multiple processes for ease of use. All scripts
work independently of each other and can be called separately. The current versions of these scripts are still in early 
development and are ***not the final version***.

## Installation
The repo can be cloned from github and the required packages with their requried version can be found in the 
"requirements.txt". In order for the application to work you need python 3.6.x or higher


To check if the installation was successful navigate to the Scripts folder and type:

`$ python3 CNV_vis.py -h` and/or `$ python3 Flag_placer.py -h`

If installation went correctly a help menu should appear in your terminal.

## Usage

There are multiple scripts for creating different tracks that can all be used independantly from each other. However the
Start_job.py is the main script that combines the tracks of the Flag_placer.py script and Softclip_graph.py script for 
each chromosome. All the scripts use the same Settings.txt file where the user can configure the settings of the 
scripts.

### Settings

- `log=True/False` Bool specifying if logfile should be made.


- `high_insert_size=-1/(any integer above 0)` Length of insert size to be classified as high or -1 if threshold should 
be calculated.
  
  
- `ultra_high_insert_size=(any integer above 0)` Length of insert size to be classified as ultra high.


- `MinPercentage_(flag name)=(any integer above 0 or 0)` float specifying thethreshold for the minimum percentage of total
  reads in region before flagged.
  
  
- `MinReadCount_(flag name)=(any integer above 0 or 0)` integer specifying the minimun amount of reads before a flag is 
  made.

### Start_job.py

The script can be called from the terminal using the following command:

`$ python3 start_job.py --bam path/to/bamfile.bam`

#### required arguments
- `--bam` or `-b` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.

#### optionalarguments

- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as Start_job.py*.
  
- `--name` or `-n` followed by `any string`. Used to give a name to the output files. *Default: "output"*
  

- `--cores` or `-c` followed by `(any integer)`. Used to specify the number of processes the script is allowed to use. 
  *Default: 1*.
  

### CNV_vis.py
The script can be called from the terminal using the following command:

`$ python3 CNV_vis.py --bam path/to/bamfile.bam --vcf_file path/to/vcffile.vcf`


#### required arguments
- `--bam` or `-b` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.
  

- `--vcf` or `-v` followed by `path/to/vcffile.vcf`. Used to specify the location of the vcf file.

#### optional arguments

- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as CNV_vis.py*.
  

- `--log` or `-l` followed by `True` or `False`. Used to specify if log file should be written in output folder. 
*Default: True*.
  
  
- `--capture_reqion` or `-cr` followed by `(any integer)`. Used to specify the amount of basepairs that should be
  included for the statistics track on the left and right of the variant calls in the vcf files. *Default: 3000*.
  
  
- `--high_insert_size` or `-hi` followed by `(any integer)`. Used to specify the threshold as what should be classified
as remarkably high insert size. *Default: 500*.
  
### Flag_placer.py
The script can be called from the terminal using the following command:

`$ python3 Flag_placer.py --bam path/to/bamfile.bam --region "chr1" `


#### required arguments
- `--bam` or `-b` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.


- `--region` or `-r` followed by `chr#:(int):(int)`. Used to specify the region the script should analyze oruse `chr#:` 
  if the script should place flags on a specific chromosome.

#### optional arguments
- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as Flag_placer.py*.
  

- `--name` or `-n` followed by `any string`. Used to give a name to the output files. *Default: "output"*

  
### softclip_graph.py

The script can be called from the terminal using the following command:

`$ python3 Flag_placer.py --bam path/to/bamfile.bam`


#### required arguments
- `--bam` or `-b` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.


- `--region` or `-r` followed by `chr#:(int):(int)`. Used to specify the region the script should analyze oruse `chr#:` 
  if the script should place flags on a specific chromosome.

#### optional arguments
- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as softclip_graph.py*.
  

- `--name` or `-n` followed by `any string`. Used to give a name to the output files. *Default: "output"*
