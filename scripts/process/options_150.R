#!/opt/conda/envs/SCAR_env/bin/Rscript

# Make functions available

library('SCAR')
library('stringr')

# Set options for processing:
 
# For fastq to bams, set bams to TRUE
# For bams to bws, set bws to TRUE
# For bamCompare to be used, set compare_bams to TRUE. 
# For bigwigCompare to be used, set compare_bgs .
# For bams to bgs to peaks, set peaks to TRUE
# For an analysis of low signal regions within identified peaks, set protection to TRUE
# For plotting with deep tools, set plot to TRUE
# For heat mapping with deep tools, set hmap to TRUE

bams <- TRUE
bws <- TRUE
compare_bams <- FALSE
compare_bgs <- TRUE
peaks <- TRUE
protection <- TRUE
plot <- TRUE
hmap <- TRUE

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
	compare=compare_bams
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
		outdir=str_c(sample_dir, 'aligned_150/'),
		min_fragment=0,
		max_fragment=150
		)
}

if (bws) {
	
	# Make tracks
	SCAR_obj <- make_tracks_opts(
	SCAR_obj,
	outdir=str_c(sample_dir,'coverage_150/'),
	compare=compare,
	comp_op='ratio',
	bin_size=1,
	genome_size = '12100000',
	center_reads = TRUE,
	skip_non_covered = FALSE,
	normalize_using='RPGC',
	extend_reads=TRUE,
	out_type='bigwig'
	)

	if (!compare) {
		SCAR_obj <- compare_bws(
		SCAR_obj,
		outdir=str_c(sample_dir,'coverage_150/'),
		comp_op='ratio',
		bin_size=1,
		skip_zeros=TRUE,
		skip_non_covered = FALSE,
		out_type='bedgraph',
		roi=NA
		)
		}	
}

if (!is.na(SCAR_obj@sample_sheet[['control_bams']]) && peaks) { 
	
	genome_dir <- pull_setting(SCAR_obj, 'genome_dir')
	# Make bedgraphs
	SCAR_obj <- make_bgs(
		SCAR_obj,
		outdir=(str_c(sample_dir, 'bedgraphs_150/')),
		pair_lr = TRUE,
		frag_size = FALSE,
		chrom_file =  'genome/sacCer_chr_sorted.txt'
		)
	
	# Call peaks
	SCAR_obj <- call_peaks_SEACR(
		SCAR_obj,
		outdir=str_c(sample_dir, 'peaks_150/'),
		norm=FALSE,
		stringent=TRUE,
		sep=""
		)
}	
	


if (protection) {
SCAR_obj <- pf_analysis(
	SCAR_obj,
	outdir=str_c(sample_dir,'peaks_150/'),
	in_bam=NA,
	in_bg=str_c(sample_dir, 'coverage_150/', *_control.cov),
	bed_file=str_c(sample_dir, 'peaks_150/', *.bed)
	)
}

SCAR_obj <- make_matrix(
	SCAR_obj,
	outdir=str_c(sample_dir,'plots/'),
	primary='reference-point',
	s_n_c='c',
	in_str=NA,
	regions=NA
	)

if (plot) {
SCAR_obj <- plot_profile(
	SCAR_obj,
	outdir=str_c(sample_dir,'plots/'),
	matrix=NA,
	plot_name='150_plot.png',
	plot_opts=NA,
	clust=3
	)
}

if (hmap) {
SCAR_obj <- plot_profile(
	SCAR_obj,
	outdir=str_c(sample_dir,'plots/'),
	matrix=NA,
	plot_name='150_hmap.png',
	plot_opts=NA,
	clust=3
	)
}