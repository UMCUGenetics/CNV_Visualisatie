# CNV_Visualisation

## Description

The CNV_visualisatie tool is a python script used to generate a BED file containing various
statistics about reads from a bam file around regions specified in a vcf file.
This version of the script is still in early development and is ***not the final version***.

## Installation
The repo can be cloned from github and the required packages can be installed by opening your terminal
 and navigating to the folder containing "requirements.txt". Use the following command in your terminal to install the
required packages by typing:

`$ pip install -r requirements.txt`

To check if the installation was successful navigate to the Scripts folder and type:

`$ python3 CNV_vis.py -h`

If installation went correctly a help menu should appear in your terminal.

## Usage
The script can be called from the terminal using the following command:

`$ python3 CNV_vis.py --bamfile path/to/bamfile.bam --vcf_file path/to/vcffile.vcf`


#### required arguments:
- `--bamfile` or `-bf` followed by `path/to/bamfile.bam`. Used to specify the location of the bam file.
  

- `--vcf_file` or `-vcf` followed by `path/to/vcffile.vcf`. Used to specify the location of the vcf file.

#### optional arguments

- `--output` or `-o` followed by `path/to/output_folder`. Used to specify the folder where the output files should be stored.
*Default: same folder as CNV_vis.py*.
  

- `--logfile` or `-l` followed by `True` or `False`. Used to specify if log file should be written in output folder. 
*Default: True*.
  
  
- `--capture_reqion` or `-cr` followed by `(any integer)`. Used to specify the amount of basepairs that should be
  included for the statistics track on the left and right of the variant calls in the vcf files. *Default: 3000*.
  
  
- `--high_insert_size` or `-hi` followed by `(any integer)`. Used to specify the threshold as what should be classified
as remarkably high insert size. *Default: 500*.