#!/opt/conda/envs/SCAR_env/bin/Rscript

# Make functions available

library('SCAR')
library('stringr')

# Set options for processing; commonly changed options are listed here, but more are available within the functions below

bams <- FALSE #process FASTQ to BAMs with bowtie2
sf_low <- 0 #lower size limit for small fragment (SF) alignments
sf_high <- 120 #upper size limit for SF
lf_low <- 140 #lower size limit for large fragment (LF) alignments
lf_high <- 180 #upper size limit for LF
ff_low <- 0 #lower size limit for full fragment (FF) alignments
ff_high <- 200 #upper size limit for FF
af_low <- 0 #lower size limit for all fragment (AF) alignments; useful for aligning longer fragments detected in ChIP-seq
af_high <- 500 #upper size limit for AF
mode <- 'end-to-end' #'end-to-end' or 'local' alignment modes for bowtie2
mapq_val <- 30 #set the value for samtools to filter reads by alignment quality (MAPQ) score
comp_op <- 'ratio' #the mathematical operation to use with bamCompare; can be 'ratio', 'log2', 'add', 'subtract'. See deeptools2 documentation for more info
bws <- FALSE #process BAM files to bigwigs
cent <- TRUE #TRUE means use fragment centers to calculate coverage; FALSE means use fragment ends
compare_cov <- FALSE #make coverage comparison calculations with bamCompare using the comp_op function defined above
peaks <- FALSE #perform peak calling on BAM files
PCF <- 'macs' #peak calling function to use; MACS is the current standard. See GitHub macs3-project/MACS for more info. SEACR is also available 
broad_flag <- FALSE #use the built-in expanded broad peak calling limit for MACS3
stringent <- TRUE #stringent peak calling cutoff designation for SEACR
protection <- FALSE #[experimental] calculate central low levels of protection
plot <- TRUE #use deeptools2 plotProfile to plot average signal based on input BED regions
hmap <- TRUE #use deeptools2 plotHeatmap to plot signal input BED regions
map_type <- "'heatmap and colorbar'" #type of information to be given on heatmap
#plot_comp <- TRUE #[deprecated] plot comparison traces or individual coverage. Superceded by s_n_c
region_file <- TRUE #use a BED file given at command line
sort_op <- 'descend' #how to sort the data when presenting it in the heatmap
sort_type <- 'max' #type of sorting calculation to use
matrix_b <- 500 #lower limit for computeMatrix, used to map coverage onto input BED regions
matrix_a <- 500 #upper limit for computeMatrix
s_n_c <- 'c' #type of bigwig file to use for computeMatrix. 'c' uses the built in coverage comparisons. Alternatively, a path to a bigwig file can be given
cov_out <- 'bigwig' #output format for bamCoverage - 'bigwig' or 'bedgraph'
comp_out <- 'bigwig' #output format for bamCompare - 'bigwig' or 'bedgraph'. Should be the same as cov_out
norm_type <- 'RPGC' #the type of normalization to be used when calculating coverage with bamCompare. 1x normalization to genome coverage (RPGC) is recommended to deal with differences in sequencing depth
clust_val <- NA #number of clusters to be fed to deeptools2 implementation of k-means. Produces clustered heatmaps and profiles
temp_dir <- 'temp/' #a temporary directory for large files to be stored during processing by deeptools2

# Pull input variables from Singularity
env_vars <- Sys.getenv(c("sample_dir", "samples_file", "plot_regions", "job_name", "frag_set"), names=TRUE)

# Make sure env_vars can be accessed/assigned
plot_regions <- env_vars[['plot_regions']]
samples_file <- env_vars[['samples_file']]
sample_dir <- env_vars[['sample_dir']]
job_name <- env_vars[['job_name']]
run_type <- env_vars[['frag_set']]
rel_dir <- getwd()

#Set the directory tag based on input fragment type
dir_tag <- str_c('_', if (run_type == 'large') {'lf'} else if (run_type == 'small') {'sf'} else if (run_type == 'full') {'ff'} else if (run_type == 'all') {'af'})

#Set the read tag denoting if fragment ends or centers are being used
if (!is.na(cent) && (cent == TRUE)) {read_tag <- str_c('_cent')}
if (is.na(cent)) {read_tag <- str_c('_bin')}

# Create sample sheet from input sample file, can contain just sample FASTQs or sample and control
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

#Record input BED file if given for use with deeptools2

if (!(region_file)) {
	region_file <- NA
	}

if (!is.na(region_file)) {
	region_file <- str_c(rel_dir, '/', plot_regions)
	}

#Make sure BAM names are added to SCAR object if not being processed
if (!bams) {
	SCAR_obj <- add_bams(
		SCAR_obj,
		alignment_dir=str_c(sample_dir, '/aligned', dir_tag, '/')
		)
	}

#Make sure bigWig names are added to SCAR object if not being processed
if (!bws) {
	if (plot || hmap) { #if plots or heatmaps are requested

		if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE) #make sure temp dir is present for deeptools2
  		Sys.setenv(TMPDIR=temp_dir)
		
		SCAR_obj <- set_settings(SCAR_obj, compare = FALSE) #make sure comparison isn't triggered before adding names to SCAR object
		SCAR_obj <- set_settings(SCAR_obj,
			alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag, '/')) #add directory to SCAR object containing bigWigs for processing

		cov_dir <- pull_setting(SCAR_obj, 'alignment_dir') #make bigWig directory current processing directory

		SCAR_obj <- add_cov( #add current processing directory to SCAR object as coverage directory
			SCAR_obj,
			alignment_dir=cov_dir,
			comp_op=comp_op)
		}
}


if ((!compare_cov) && (plot || hmap || (protection && comp_out == 'bedgraph'))) { 
	#if comparison is not requested, make sure bigWigCompare names are added to SCAR object
	#Also necessary if bedgraph format is requested to prevent errors in SCAR object entries

	plot_dir <- str_c(sample_dir, '/plots', dir_tag, read_tag, '_', #assemble name of plotting directory
		comp_op, '_', matrix_b, '_to_', matrix_a, '/')
	SCAR_obj <- set_settings(SCAR_obj, compare = TRUE) #indicate comparison is requested to add_cov function
	SCAR_obj <- set_settings(SCAR_obj, comp_op=comp_op) #record the operation used for comparison

	comp_op <- pull_setting(SCAR_obj, 'comp_op') #pull comp_op for specifying name of comparison directory

	if (comp_out == 'bigwig') { #set bigwig output as active
		SCAR_obj <- set_settings(SCAR_obj,
		alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag,
			'/bw_comp_', comp_op, '/'))
			}

	else if (comp_out == 'bedgraph') { #set bedgraph output as active
		SCAR_obj <- set_settings(SCAR_obj,
		alignment_dir=str_c(sample_dir, '/coverage', dir_tag, read_tag,
			'/bg_comp_', comp_op, '/'))
			}

	alignment_dir <- pull_setting(SCAR_obj, 'alignment_dir') #pull active directory from SCAR object

	SCAR_obj <- add_cov( #assign active directory to comparison slot in SCAR object
		SCAR_obj,
		alignment_dir=alignment_dir,
		comp_op=comp_op
		)
	}



if (!peaks) { #if peak calling is not requested, make sure peak directories are added to SCAR object
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

if (bams) { #FASTQ processing to binary alignments is requested
	#large, small, full, all fragment ranges are defined above and requested as command line input

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


	#apply FASTQC to sample and control FASTQs
	SCAR_obj <- fastqc( 
		SCAR_obj,
		outdir=str_c(sample_dir, '/fastqc_reports/')
		)

	#index the reference genome; bowtie2 uses Burrows-Wheeler transform 
	SCAR_obj <- bowtie2_index(
		SCAR_obj,
		outdir=str_c(rel_dir, '/genome/'),
		genome_assembly=str_c(rel_dir, '/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'),
		index_name='index/sacCer3_index'
		)


	SCAR_obj <- bowtie2_align( #align FASTQs using the indexed genome with bowtie2
		SCAR_obj,
		outdir=str_c(sample_dir, '/aligned', dir_tag, '/'),
		alignment_mode=mode,
		min_fragment=min_fragment,
		max_fragment=max_fragment,
		mapq_val=mapq_val
		)
}

if (bws) { #coverage calculation requested based on BAMs
	SCAR_obj <- set_settings(SCAR_obj, compare=FALSE) #make sure outputs are added as coverage and not comparisons
	cov_dir <- str_c(sample_dir, '/coverage', dir_tag, read_tag, '/') #build coverage directory name
	
	center_reads <- cent #record decision on whether to use fragment centers or ends for coverage calculation

	# Make tracks using input options, generally leave this alone unless you have a specific option you know needs to be changed here
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

if (compare_cov) { #comparison of coverage is requested

	SCAR_obj <- set_settings(SCAR_obj, comp_op=comp_op) #add to comparison function to SCAR obj
	comp_op <- pull_setting(SCAR_obj, 'comp_op')
	comp_dir <- if (comp_out == 'bigwig') { #make sure appropriate directory name is recorded depending on output file structure requested
		str_c(sample_dir, '/coverage', dir_tag, read_tag, 
			'/bw_comp_', comp_op, '/')
			}
		else if (comp_out == 'bedgraph') {
			str_c(sample_dir, '/coverage', dir_tag, read_tag,
				'/bg_comp_', comp_op, '/')
				}

	SCAR_obj <- compare_bws( #perform coverage comparison between sample and control bigWigs or bedgraphs
		SCAR_obj,
		outdir=comp_dir,
		comp_op=comp_op,
		bin_size=1,
		skip_zeros=FALSE,
		out_type=comp_out,
		roi=NA
		)
	}

if (!is.na(SCAR_obj@sample_sheet[['control_bams']]) && peaks) { #peak calling requested; make sure a control sample is present for peak calling

	if (PCF == 'seacr') { #if SEACR was requested as peak calling function

		genome_dir <- pull_setting(SCAR_obj, 'genome_dir') #grab genome directory for chromosome info to make bedgraph files

		#build the required bedgraph inputs for SEACR
		SCAR_obj <- make_bgs( 
			SCAR_obj,
			outdir=str_c(sample_dir, '/bedgraphs', dir_tag, '/'),
			pair_lr=TRUE,
			frag_size=FALSE,
			chrom_file='genome/sacCer_chr_sorted.txt'
			)

		# Call peaks with SEACR; also adds SEACR outputs to SCAR object
		SCAR_obj <- call_peaks_SEACR(
			SCAR_obj,
			outdir=str_c(sample_dir, '/peaks', dir_tag, '/'),
			norm=TRUE,
			stringent=stringent
			)
		
		SCAR_obj <- add_beds( #add dummy MACS outputs to satisfy SCAR object
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/macs/'),
			peak_type = 'macs',
			stringent = stringent
			)
	}

	if (PCF == 'macs') { #MACS3 peak calling requested - MACS is the standard and preferred method

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

		SCAR_obj <- add_beds( #add dummy SEACR output to satisfy SCAR object
			SCAR_obj,
			peak_dir=str_c(sample_dir, '/peaks', dir_tag, '/seacr/'),
			peak_type = 'seacr',
			stringent = stringent)	
	}
}

if (protection) { #experimental feature for identifying low coverage in the center of otherwise high signal regions

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

if (plot || hmap) { #plotting is requested using deeptools2 plotProfile or plotHeatmap

	plot_dir <- str_c(sample_dir, '/plots', dir_tag, read_tag, '/') #build plot output directory name

	id_str <- '_no_clust' #standard is to add a tag indicating no clustering

	if (!is.na(clust_val)) { #add cluster number tag if requested
		id_str <- str_c('_clust_', clust_val)
		}

	if (sort_op == 'keep') { #keep indicates no sorting and is used when making order matched heatmaps
		id_str <- str_c('_linked')
		}

	id_str <- str_c(id_str, '_', matrix_b, '_', matrix_a) #add plotting bounds to output file name
	
	if (!is.na(region_file)) { #add name identifier to output files
		id_str <- str_c(id_str, '_', job_name)
		}

	if (is.na(region_file)) { #if no region file was given, add job name as an identifier
		id_str <- str_c(id_str, '_in_sites')
		}
	}
	
	out_loc <- str_c(plot_dir, sample_id, dir_tag, '_', comp_op, id_str, read_tag, '/') #assemble plotting directory info into an output directory name
	
	SCAR_obj <- make_matrix( #map coverage onto input BED regions; will use peaks file if no BED is given and peaks were requested, otherwise will not run
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

if (plot) { #profile plotting requested
	
	out_loc <- str_c(plot_dir, sample_id, dir_tag, '_', comp_op, id_str, read_tag, '/') #redundant assignment of outdir for profile
	
	points_file <- str_c(out_loc, '/', sample_id, '_profile', dir_tag, '_', comp_op, id_str, read_tag, '.tsv') #make the name for a TSV of profile points

	SCAR_obj <- plot_profile( #create an average profile plot based on output of computeMatrix
		SCAR_obj,
		outdir=out_loc,
		matrix=NA,
		plot_name=str_c(sample_id, '_profile', dir_tag, '_', comp_op, id_str, read_tag, '.png'), #name and format of profile plot
		plot_opts=str_c(if(!is.na(clust_val)) {str_c('--kmeans', clust_val, sep=' ')}, #options based on deeptools2 docs, can include any desired bash flags
				'--plotHeight', '10',
				'--plotWidth', '10',
				'-T', "'Coverage signal in given regions'",
				'--perGroup',
				'--outFileNameData', str_c(points_file),
				sep = " "),
		clust=clust_val
		)
	}

if (hmap) { #heatmap is requested, will be output to out_loc
	SCAR_obj <- plot_heatmap(
		SCAR_obj,
		outdir=out_loc,
		matrix=NA,
		plot_name=str_c(sample_id, '_hmap', dir_tag, '_', comp_op, id_str, read_tag, '.png'),
		plot_opts=str_c(if(!is.na(clust_val)) {str_c('--kmeans', clust_val, sep=' ')}, #options based on deeptools2 docs, can include any desired bash flags
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
