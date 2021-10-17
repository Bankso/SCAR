#!/usr/bin/Rscript

#This is the R script that does the heavy lifting
#Tells all the programs to do their thing in the correct order and with the correct resource allocation

# Create sample sheet
samples <- read.delim("samples.txt", sep = "\t")

# Generate SCAR object
library("SCAR")

SCAR <- SCAR_maker(
  analysis_type = "SChEC-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)

# Perform read QC
SCAR <- fastqc(SCAR, outdir = "./fastqc_reports")

# Generate Bowtie2 index and align reads
SCAR <- bowtie2_index(
  SCAR, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  index_name = "sacCer3_index"
)

SCAR <- bowtie2_align(SCAR, outdir = "./aligned", min_fragment = 20, max_fragment = 200)

# Make tracks
SCAR <- make_bigwigs(
  SCAR, outdir = "./bigwigs", bin_size = 1,
  normalize_using = "CPM", extend_reads = TRUE
)

# Call peaks
SCAR <- call_peaks(SCAR, outdir = './peaks', genome_size = 12100000)
