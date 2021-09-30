#!/usr/bin/Rscript

#This is the R script that does the heavy lifting
#Tells all the programs to do their thing in the correct order and with the correct resource allocation

# Create sample sheet
samples <- read.delim("samples.txt", sep = "\t")

# Generate ZentTools object
library("ZentTools")

zent <- zent_tools(
  analysis_type = "ChEC-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)

# Perform read QC
zent <- fastqc(zent, outdir = "./fastqc_reports")

# Generate Bowtie2 index and align reads
zent <- bowtie2_index(
  zent, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  index_name = "sacCer3_index"
)

zent <- bowtie2_align(zent, outdir = "./aligned", min_fragment = 3, max_fragment = 100)

# Make tracks
zent <- make_bigwigs(
  zent, outdir = "./bigwigs", bin_size = 1,
  normalize_using = "CPM", extend_reads = TRUE
)

# Call peaks
zent <- call_peaks(zent, outdir = './peaks', genome_size = 12100000)
