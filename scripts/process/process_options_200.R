#!/opt/conda/envs/SCAR_env/bin/Rscript

# Make functions available

library('SCAR')
library('stringr')

# Set options for processing:
 
# For fastq to bams, set bams to TRUE
# For bams to bws, set bws to TRUE
# For bams to bgs to peaks, set peaks to TRUE
# For an analysis of low signal regions within identified peaks, set protection to TRUE

bams <- TRUE
bws <- TRUE
peaks <- TRUE
protection <- TRUE

# Pull input variables from Singularity
env_vars <- Sys.getenv(c("sample_dir", "samples_file"), names=TRUE)

# Make sure env_vars can be accessed/assigned
samples_file <- env_vars[['samples_file']]
sample_dir <- env_vars[['sample_dir']]
rel_dir <- getwd()

# Create sample sheet from input sample file
samples <- read.delim(samples_file, sep='\t')

# Create SCAR_obj, holds settings and file paths
SCAR_obj <- SCAR_maker(
	analysis_type='SChEC-seq', 
	sample_sheet=samples,
	paired=TRUE, 
	ncores=8, 
	compare=TRUE
	)

if (bams) {
	
	# Perform read QC
	SCAR_obj <- fastqc(
		SCAR_obj, 
		outdir=str_c(sample_dir, 'fastqc_reports/')
		)

	# Create the bowtie2 genome index from sacCer3
	SCAR_obj <- bowtie2_index(
		SCAR_obj,
		outdir=str_c(rel_dir, '/genome/'),
		genome_assembly=str_c(rel_dir, '/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'),
		index_name='index/sacCer3_index'
		)


	SCAR_obj <- bowtie2_align(
		SCAR_obj,
		outdir=str_c(sample_dir, 'aligned_200/'),
		min_fragment=0,
		max_fragment=150
		)
}

else {
	set_settings(SCAR_obj, analysis_dir = str_c(sample_dir, 'aligned/'))
	
	add_bams(SCAR_obj, analysis_dir = analysis_dir)
if (bws) {
	
	# Make tracks
	SCAR_obj <- make_tracks(
	SCAR_obj,
	outdir=str_c(sample_dir,'coverage_200/'),
	comp_op='ratio',
	bin_size=1,
	genome_size = '12100000',
	center_reads = TRUE,
	skip_non_covered = FALSE,
	normalize_using='RPGC',
	extend_reads=TRUE
	)	
}

if (!is.na(SCAR_obj@sample_sheet[['control_bams']]) && peaks) { 
	
	genome_dir <- pull_setting(SCAR_obj, 'genome_dir')
	# Make bedgraphs
	SCAR_obj <- make_bgs(
		SCAR_obj,
		outdir=(str_c(sample_dir, 'bedgraphs_200/')),
		pair_lr = TRUE,
		frag_size = FALSE,
		chrom_file =  'genome/sacCer_chr_sorted.txt'
		)
	
	# Call peaks
	SCAR_obj <- call_peaks_SEACR(
		SCAR_obj,
		outdir=str_c(sample_dir, 'peaks_200/'),
		norm=FALSE,
		stringent=TRUE,
		sep=""
		)
}	
	


if (protection) {
SCAR_obj <- pf_analysis(
	SCAR_obj,
	outdir=str_c(sample_dir,'peaks_200/'),
	bed_file = str_c(sample_dir, 'peaks_200/', *.bed)
	)
}