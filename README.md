# Mock Home
### A repository for a directory framework that fits the SCAR (SpLiT-ChEC Analysis with R) bioinformatics pipeline
See the READMEs in [SCAR](https://github.com/Bankso/SCAR) and [SEAPE](https://github.com/Bankso/SEAPE) for more information.

#### How to prepare *mock_home*
After downloading, place *mock_home/* in your home directory (or somewhere you have permissions).

Loading in raw sequencing reads:  
1) Place your sample FASTQ(s) into *home_dir/samples/*  
2) Place your control FASTQ(s) into *home_dir/controls/*   
3) Edit *home_dir/samples/samples.txt* to fit your sample and control names for processing, as well as file paths for single or paired-end reads for each sample (compressed or uncompressed)

Verifying start script and settings:
1) Open *mock_home/scripts/start/run.sh* in a text editor  
2) Change the batch submission script to fit your submission architecture/details  
3) Save the start script to *mock_home/scripts/start/* - if the name was changed, use that as the *$path/to/run* below
4) Modify any settings you want to change in options.R and save
  
### Starting a processing run
To initiate a run from the command line, use the following call format:
```
sbatch $path/to/run.sh $path_to_new_home $path/to/options.R $name_for_output $path_to_sample_files $path_to_BED_file $fragment_range
```
Where $fragment_range can be one of four settings, each of which signals a different range of values indicating fragment sizes to be processed during alignment to a reference genome: all, full, small, and large. Default ranges are 0-500bp, 0-200bp, 0-120bp, and 140-180bp, respectively. This input is used in conjunction with a set of output type identifiers (noted below) to create identifiable directories from variably processed data, simplifying downstream analysis and process optimization.

Example command:
```
sbatch scripts/start/run.sh scripts/process/options.R ChIP_FLAG Abf1/rep1 BED_files/plot_regions.bed full
```
1) Navigate to the folder *mock_home/*
2) From the command line, enter the batch submission call for your HPC with the desired input parameters to initiate processing
3) Genome database files and a Singularity image file will be downloaded to the mock_home directory before processing begins
4) Files will be output to *mock_home/samples/outdir* for each function - output files will have names roots corresponding to sample or control names in *samples.txt*, as well as multiple identifiers that simplify directory organization and identification of settings associated with outputs. (Detailed list of identifiers in progress)

### Output files and directory structure

**All files can be found in *mock_home/samples/***  
*fastqc_reports/* - fastq quality reports from FASTQC in HTML format, easily viewable in a standard browser.  
*aligned/* - SAM, BAM, BAM indexes from bowtie2 and samtools  
*coverage/* - bigWig or bedgraph coverage files from deepTools   
*peaks/* - peak finding outputs from MACS3  
*plots/* - plots, matrices, heatmaps, and sorted region files from plotting functions  
