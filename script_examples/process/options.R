#!/opt/conda/envs/SCAR_env/bin/Rscript

# Make functions available

library('SCAR')
library('stringr')

# Set options for processing:

bams <- FALSE		
sf_low <- 0
sf_high <- 120
lf_low <- 140
lf_high <- 180
ff_low <- 0
ff_high <- 200
af_low <- 0
af_high <- 500
mode <- 'end-to-end'
mapq_val <- 30
comp_op <- 'ratio'
bws <- FALSE
cent <- TRUE
compare_cov <- FALSE
peaks <- FALSE
PCF <- 'macs'
broad_flag <- FALSE
stringent <- TRUE
protection <- FALSE
plot <- TRUE
hmap <- TRUE
map_type <- "'heatmap and colorbar'"
plot_comp <- TRUE
region_file <- TRUE
sort_op <- 'descend'
sort_type <- 'max'
matrix_b <- 500
matrix_a <- 500
s_n_c <- 'c'
cov_out <- 'bigwig'
comp_out <- 'bigwig'
norm_type <- 'RPGC'
clust_val <- NA
temp_dir <- 'temp/'

# Pull input variables from Singularity
env_vars <- Sys.getenv(c("sample_dir", "samples_file", "plot_regions", "job_name", "frag_set"), names=TRUE)

# Make sure env_vars can be accessed/assigned
plot_regions <- env_vars[['plot_regions']]
samples_file <- env_vars[['samples_file']]
sample_dir <- env_vars[['sample_dir']]
job_name <- env_vars[['job_name']]
run_type <- env_vars[['frag_set']]
rel_dir <- getwd()

dir_tag <- str_c('_', if (run_type == 'large') {'lf'} else if (run_type == 'small') {'sf'} else if (run_type == 'full') {'ff'} else if (run_type == 'all') {'af'})

if (!is.na(cent) && (cent == TRUE)) {read_tag <- str_c('_cent')}
if (is.na(cent)) {read_tag <- str_c('_bin')}

# Create sample sheet from input sample file
samples <- read.delim(samples_file, sep='\t')

# Create SCAR_obj, holds settings and file paths
SCAR_obj <- SCAR_maker(
	analysis_type='SChEC-seq',
	sample_sheet=samples,
	paired=TRUE,
	ncores=8,
	compare=FALSE,
	comp_op=comp_op
	)

sample_id <- SCAR_obj@sample_sheet[['sample_name']]

# Logic gates to decide how to name files and what to process

if (!(region_file)) {
	region_file <- NA
	}

if (!is.na(region_file)) {
	region_file <- str_c(rel_dir, '/', plot_regions)
	}

if (!bams) {
	SCAR_obj <- add_bams(
		SCAR_obj,
		alignment_dir=str_c(sample_dir, '/aligned', dir_tag, '/')
		)
	}

if (!bws) {
	if (plot || hmap) {

		if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  		Sys.setenv(TMPDIR=temp_dir)
		
		SCAR_obj <- set_settings(SCAR_obj, compare = FALSE)
		SCAR_obj <- set_settings(SCAR_obj,
			alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag, '/'))

		cov_dir <- pull_setting(SCAR_obj, 'alignment_dir')

		SCAR_obj <- add_cov(
			SCAR_obj,
			alignment_dir=cov_dir,
			comp_op=comp_op)
		}
}

if ((!compare_cov) && (plot || hmap || (protection && comp_out == 'bedgraph'))) {

	plot_dir <- str_c(sample_dir, '/plots', dir_tag, read_tag, '_',
		comp_op, '_', matrix_b, '_to_', matrix_a, '/')
	SCAR_obj <- set_settings(SCAR_obj, compare = TRUE)
	SCAR_obj <- set_settings(SCAR_obj, comp_op=comp_op)

	comp_op <- pull_setting(SCAR_obj, 'comp_op')

	if (comp_out == 'bigwig') {
		SCAR_obj <- set_settings(SCAR_obj,
		alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag,
			'/bw_comp_', comp_op, '/'))
			}

	else if (comp_out == 'bedgraph') {
		SCAR_obj <- set_settings(SCAR_obj,
		alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag,
			'/bg_comp_', comp_op, '/'))
			}

	alignment_dir <- pull_setting(SCAR_obj, 'alignment_dir')

	SCAR_obj <- add_cov(
		SCAR_obj,
		alignment_dir=alignment_dir,
		comp_op=comp_op
		)
	}



if (!peaks) {
	SCAR_obj <- add_beds(
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/seacr/'),
			peak_type = 'seacr',
			stringent = stringent)

	SCAR_obj <- add_beds(
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/macs/'),
			peak_type = 'macs',
			stringent = stringent)
	}

if (bams) {

	if (run_type == 'large') {
		min_fragment <- lf_low
		max_fragment <- lf_high
	}
	
	if (run_type == 'small') {
		min_fragment <- sf_low
		max_fragment <- sf_high
	}

	if (run_type == 'full') {
		min_fragment <- ff_low
		max_fragment <- ff_high
	}

	if (run_type == 'all') {
		min_fragment <- af_low
		max_fragment <- af_high
	}


	# Perform read QC
	SCAR_obj <- fastqc(
		SCAR_obj,
		outdir=str_c(sample_dir, '/fastqc_reports/')
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
		outdir=str_c(sample_dir, '/aligned', dir_tag, '/'),
		alignment_mode=mode,
		min_fragment=min_fragment,
		max_fragment=max_fragment,
		mapq_val=mapq_val
		)
}

if (bws) {
	SCAR_obj <- set_settings(SCAR_obj, compare=FALSE)
	cov_dir <- str_c(sample_dir, '/coverage', dir_tag, read_tag, '/')
	
	center_reads <- cent

	# Make tracks
	SCAR_obj <- make_tracks_opts(
	SCAR_obj,
	outdir=cov_dir,
	compare=FALSE,
	MNase=FALSE,
	comp_op=comp_op,
	bin_size=1,
	genome_size=12100000,
	center_reads=center_reads,
	min_fragment=NA,
	max_fragment=NA,
	normalize_using=norm_type,
	extend_reads=TRUE,
	out_type=cov_out
	)
	}

if (compare_cov) {

	SCAR_obj <- set_settings(SCAR_obj, comp_op=comp_op)
	comp_op <- pull_setting(SCAR_obj, 'comp_op')
	comp_dir <- if (comp_out == 'bigwig') {
		str_c(sample_dir, '/coverage', dir_tag, read_tag, 
			'/bw_comp_', comp_op, '/')
			}
		else if (comp_out == 'bedgraph') {
			str_c(sample_dir, '/coverage', dir_tag, read_tag,
				'/bg_comp_', comp_op, '/')
				}

	SCAR_obj <- compare_bws(
		SCAR_obj,
		outdir=comp_dir,
		comp_op=comp_op,
		bin_size=1,
		skip_zeros=FALSE,
		out_type=comp_out,
		roi=NA
		)
	}

if (!is.na(SCAR_obj@sample_sheet[['control_bams']]) && peaks) {

	if (PCF == 'seacr') {

		genome_dir <- pull_setting(SCAR_obj, 'genome_dir')

		# Make bedgraphs
		SCAR_obj <- make_bgs(
			SCAR_obj,
			outdir=str_c(sample_dir, '/bedgraphs', dir_tag, '/'),
			pair_lr=TRUE,
			frag_size=FALSE,
			chrom_file='genome/sacCer_chr_sorted.txt'
			)

		# Call peaks
		SCAR_obj <- call_peaks_SEACR(
			SCAR_obj,
			outdir=str_c(sample_dir, '/peaks', dir_tag, '/'),
			norm=TRUE,
			stringent=stringent
			)
		
		SCAR_obj <- add_beds(
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/macs/'),
			peak_type = 'macs',
			stringent = stringent
			)
	}

	if (PCF == 'macs') {

		SCAR_obj <- call_peaks_macs(
			SCAR_obj,
			outdir = str_c(sample_dir, '/peaks', dir_tag, '/'),
			process = "callpeak",
			in_form = "BAMPE",
			genome_size = 12100000,
			tag_size = NA,
			width = NA,
			mfold = NA,
			q_val = NA,
			p_val = NA,
			broad = broad_flag
			)

		SCAR_obj <- add_beds(
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/seacr/'),
			peak_type = 'seacr',
			stringent = stringent)	
	}
}

if (protection) {

# Put together inputs for low signal analysis
# If you want to use an input type, set it to the named variable in the
# 'pf_analysis()' function below this input section. Select only one input type.


# Change to desired bam file path if required

in_bam <- SCAR_obj@sample_sheet[['sample_bams']]

in_bam <- in_bam


# Change to desired coverage file, bedgraph format required

in_bg <- SCAR_obj@sample_sheet[['compared_cov']]

in_bg <- NA


# Change to desired bed file path if required

peak_bed <- str_c(PCF, '_peaks')

bed_file <- SCAR_obj@sample_sheet[[peak_bed]]

bed_file <- bed_file


SCAR_obj <- pf_analysis(
	SCAR_obj,
	outdir=str_c(sample_dir, '/peaks', dir_tag, '/'),
	in_bam=in_bam,
	in_bg=in_bg,
	thresh=2,
	bed_file=bed_file
	)
}

if (plot || hmap) {

	plot_dir <- str_c(sample_dir, '/plots', dir_tag, read_tag, '/')

	id_str <- '_no_clust'

	if (!is.na(clust_val)) {
		id_str <- str_c('_clust_', clust_val)
		}

	if (sort_op == 'keep') {
		id_str <- str_c('_linked')
		}

	id_str <- str_c(id_str, '_', matrix_b, '_', matrix_a)
	
	if (!is.na(region_file)) {
		id_str <- str_c(id_str, '_', job_name)
		}

	if (is.na(region_file)) {
		id_str <- str_c(id_str, '_in_sites')
		}
	}
	
	out_loc <- str_c(plot_dir, sample_id, dir_tag, '_', comp_op, id_str, read_tag, '/')
	
	SCAR_obj <- make_matrix(
		SCAR_obj,
		outdir=out_loc,
		primary='reference-point',
		s_n_c=s_n_c,
		in_str=str_c('--referencePoint', 'center',
			'-b', matrix_b,
			'-a', matrix_a,
			'-bs', 1,
			sep = ' '),
		regions=region_file,
		peak_type=PCF
		)

if (plot) {
	
	out_loc <- str_c(plot_dir, sample_id, dir_tag, '_', comp_op, id_str, read_tag, '/')
	
	points_file <- str_c(out_loc, '/', sample_id, '_profile', dir_tag, '_', comp_op, id_str, read_tag, '.tsv')

	SCAR_obj <- plot_profile(
		SCAR_obj,
		outdir=out_loc,
		matrix=NA,
		plot_name=str_c(sample_id, '_profile', dir_tag, '_', comp_op, id_str, read_tag, '.png'),
		plot_opts=str_c(if(!is.na(clust_val)) {str_c('--kmeans', clust_val, sep=' ')},
				'--plotHeight', '10',
				'--plotWidth', '10',
				'-T', "'Coverage signal in given regions'",
				'--perGroup',
				'--outFileNameData', str_c(points_file),
				sep = " "),
		clust=clust_val
		)
	}

if (hmap) {
	SCAR_obj <- plot_heatmap(
		SCAR_obj,
		outdir=out_loc,
		matrix=NA,
		plot_name=str_c(sample_id, '_hmap', dir_tag, '_', comp_op, id_str, read_tag, '.png'),
		plot_opts=str_c(if(!is.na(clust_val)) {str_c('--kmeans', clust_val, sep=' ')},
				'-T', "'Coverage signal in given regions'",
				'--colorMap', 'Greens',
				'--boxAroundHeatmaps', 'no',
				'--sortRegions', sort_op,
				'--sortUsing', sort_type,
				'--whatToShow', map_type,
				sep = ' '),
			clust=clust_val
		)
	}
