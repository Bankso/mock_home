# Mock Home
## A repository for the file framework used in data processing via SCAR (SpLiT-ChEC Analysis with R)  
See the README in *'Bankso/SCAR'* for more information.  
Not intended to be used independently of SCAR processing.

### How to use
After downloading, place *mock_home/* in your home directory (or somewhere you have permissions).

Input file setup:  
1) place your sample fastqs into *home_dir/samples/*  
2) place your control fastqs into *home_dir/controls/*   
3) edit *home_dir/samples/samples.txt* to fit your sample, control, and fastq names 
  
Start script setup:  
1) Open *mock_home/scripts/start/q_run.sh* in a text editor  
2) Change the batch submission script to fit your submission architecture/details  
3) Save the start script to *mock_home/scripts/start/* - if the name was changed, use that as the *run_script* name below 
  
Start processing run:  
command format: *sbatch* *run script* *options script* *run_name*  
1) Navigate to the folder *mock_home/*
2) From the command line, start with the batch submission call for your HPC (*i.e.*, *sbatch*) 
2) paste the following after the batch call (edit to fit your details): *scripts/start/q_run.sh scripts/process/options_#.R new_run_1*  
	 
Example command: *sbatch scripts/start/q_run.sh scripts/process/options_200.R Abf1_SC_5s*  

Files will be output to *mock_home/samples/outdir* for each function - output files will have names corresponding to sample and control names in *samples.txt*  

### Output files and directory structure

**All files can be found in *mock_home/samples/***  
*fastqc/* - fastq quality reports from FASTQC  
*aligned/* - SAM, BAM, BAM indexes from bowtie2 and samtools  
*coverage/* - bigWig or bedgraph coverage files from deepTools  
*bedgraphs/* - bedgraph coverage files from bedTools - used as input for peak finding  
*peaks/* - peak bed and overlap bed region files from SEACR and low signal finder  
*plots/* - plots, matrices, heatmaps, and sorted region files from plotting functions  
