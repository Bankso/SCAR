#!/usr/bin/Rscript

#This is the R script that does the heavy lifting
#Tells all the programs to do their thing in the correct order and with the correct resource allocation

# Create sample sheet
samples <- read.delim("samples.txt", sep = "\t")

# Generate SCAR_obj
library("SCAR")

SCAR_obj <- SCAR_maker(
  analysis_type = "SChEC-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)

# Perform read QC
SCAR_obj <- fastqc(SCAR_obj, outdir = "./fastqc_reports")

# Generate Bowtie2 index and align reads
SCAR_obj <- bowtie2_index(
  SCAR_obj, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  index_name = "sacCer3_index"
)

SCAR_obj <- bowtie2_align(SCAR_obj, outdir = "./aligned", min_fragment = 20, max_fragment = 200)

# Make tracks
SCAR_obj <- make_bigwigs(
  SCAR_obj, outdir = "./bigwigs", bin_size = 1,
  normalize_using = "CPM", extend_reads = TRUE
)

# Call peaks
SCAR_obj <- call_peaks_SEACR(SCAR_obj, outdir = './peaks')
