# CNV_Visualisation

## Description

The CNV_visualisation project is a project about visually detecting and analyzing CNV calls. It contains multiple python
scripts that retrieve additional information from bam files. The CNV_vis.py script creates a BED formatted file that can
be loaded in IGV and gives a summary of the read information around CNV calls made in a vcf file. The flag_placer.py is
a script creates a BED formatted file that highlights potential regions of interest in CNV calling. Both scripts work
independently of each other and can be called separately.
The current versions of these scripts are still in early development and are ***not the final version***.

## Installation
The repo can be cloned from github and the required packages with their requried version can be found in the "requirements.txt".


To check if the installation was successful navigate to the Scripts folder and type:

`$ python3 CNV_vis.py -h` and/or `$ python3 Flag_placer.py -h`

If installation went correctly a help menu should appear in your terminal.

## Usage
###CNV_vis.py
The script can be called from the terminal using the following command:

`$ python3 CNV_vis.py --bam path/to/bamfile.bam --vcf_file path/to/vcffile.vcf`


#### required arguments:
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

`$ python3 Flag_placer.py --bam path/to/bamfile.bam`


#### required arguments
- `--bam` or `-b` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.

#### optional arguments
- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as Flag_placer.py*.
  

- `--log` or `-l` followed by `True` or `False`. Used to specify if log file should be written in output folder. 
*Default: True*.
  

- `--threshold` or `-t` followed by `(any integer)`. Used to specify the minimum number of reads with a relevant feature
before it is flagged. *Default: 0*.
  

- `--minpercentage` or `-mp` followed by `(any float)`. Used to specify the minimum percentage of reads in a region of
interest before it is marked as a region of interest. *Default: 0*
  

- `--region` or `-r` followed by `chr#:(int):(int)`. Used to specify the region the script should analyze. Use `all` if 
the script should place flags on the whole genome or use `chr#:0-0` if the script should place flags on a specific 
  chromosome. *Default: "all"*
  
- `--high_insert_size` or `-hi` followed by `(any integer)`. Used to specify the threshold as what should be classified
as remarkably high insert size. *Default: 500*.
  

