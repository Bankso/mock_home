# Mock Home
## A repository for the file framework used in data processing via SCAR (SpLiT-ChEC Analysis with R)

See the README in Bankso/SCAR for more information.  
Not intended to be used independently of SCAR processing.

After downloading, place the mock_home file in your home directory (or somewhere you have permissions).

Input file setup:  
1) place your sample fastqs into samples/  
2) place your control fastqs into controls/   
3) edit 'samples.txt' in samples/ to fit your sample, control, and fastq names.  
  
Start script setup:  
1) Open mock_home/scripts/start/run_start.sh in a text editor  
2) Change the batch submission script to fit your submission architecture/details  
3) Save the start script to mock_home/scripts/start/  
4) Navigate to the mock_home directory  
  
Start processing run:  
1) From the command line, while inside mock_home/, start with the batch submission call for your HPC  
	i.e., sbatch  
2) paste the following after the batch call (edit to fit your details):  
	 scripts/start/run_start_2.sh scripts/process/process_options_#.R new_run_1  
  
Files should be dumped to samples/ - outputs should have names corresponding to those in process_options