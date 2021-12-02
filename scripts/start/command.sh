#!/usr/bin/bash

#SBATCH --job-name=SCAR	    		 	### Job Name
#SBATCH --output=scripts/logs/out.txt           ### File in which to store job output DO NOT EDIT
#SBATCH --error=scripts/logs/err.txt            ### File in which to store job error messages DO NOT EDIT
#SBATCH --time=0-01:00:00           		### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   		### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         		### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8	    		### Number of CPU's used to complete submitted tasks within each requested node
#SBATCH --mem=12GB		    		### Amount of preallocated memory for submitted job
#SBATCH --account=mcknightlab       		### Account used for job submission

## Usage: sbatch command_start.sh 1 2 3
## Requires two inputs at the command line:
	## For start bash script:
		## First input - home directory for data processing
		## Second input - function call at command line (i.e., bash, bamCompare, etc.) with args
		## Third input - what to name output(i.e., 'new_file', etc.)
home_dir=$1
input=$2
out=$3

cd $home_dir 

module load singularity
singularity exec --env input='$input',out='$out' -eCB $home_dir:/opt/conda/process -H /opt/conda/process *.sif $input
mv scripts/logs/out.txt $home_dir/$out

# Clean-up from processing
rm -r ./temp/
rm ./SEACR_1.3.*

exit
