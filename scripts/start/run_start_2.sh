#!/usr/bin/bash

#SBATCH --job-name=SCAR	    		 	### Job Name
#SBATCH --output=scripts/logs/run.out           ### File in which to store job output DO NOT EDIT
#SBATCH --error=scripts/logs/run.err            ### File in which to store job error messages DO NOT EDIT
#SBATCH --time=0-01:00:00           		### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   		### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         		### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8	    		### Number of CPU's used to complete submitted tasks within each requested node
#SBATCH --mem=12GB		    		### Amount of preallocated memory for submitted job
#SBATCH --account=your_account_name       		### Account used for job submission

## Be sure to edit the above to fit your HPC job submisssion framework
## For example directory structure, see the README at https://github.com/Bankso/SCAR
## Usage: sbatch all_start.sh 1 2
## Requires two inputs at the command line:
	## For start bash script:
		## First input - Path to R script for processing
		## Second input - Name of job (i.e., Gal4_chec or something like that) for sbatch log files
		
home_dir=$PWD
process_script=$1
job_name=$2
sample_dir=$home_dir/samples
samples_file=$sample_dir/samples.txt

cd $home_dir

mkdir -p genome && cd genome

wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gff3.gz

gunzip Saccharomyces_cerevisiae.R64-1-1.100.gff3.gz

cd $home_dir

module load singularity

singularity cache clean -f

if [ -e ./scar* ]; then
	rm scar*
	singularity pull --arch amd64 library://banksorion/default/scar_software:latest_v1.*
else
	singularity pull --arch amd64 library://banksorion/default/scar_software:latest_v1.*
fi

# Start singularity-based processing run 
singularity exec --env sample_dir=$sample_dir,samples_file=$samples_file -eCB $home_dir:/opt/conda/process -H /opt/conda/process *.sif Rscript $process_script

# Move sbatch logs to sample dir for storage
mkdir -p $sample_dir/logs/ && mv scripts/logs/run.out $sample_dir/logs/$job_name.out && mv scripts/logs/run.err $sample_dir/logs/$job_name.err

# Clean-up from processing
rm -r ./temp/
rm ./SEACR_1.3.*

exit