% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setting_handler.R
\name{set_settings}
\alias{set_settings}
\title{Set Settings.}
\usage{
set_settings(
  SCAR_obj,
  analysis_type = NA,
  paired = NA,
  ncores = NA,
  genome_dir = NA,
  genome_index = NA,
  genome_assembly = NA,
  alignment_dir = NA,
  peak_dir = NA,
  compare = NA,
  plot_dir = NA,
  comp_op = NA
)
}
\arguments{
\item{SCAR_obj}{SCAR object, holds data for use in all the functions.}

\item{analysis_type}{Type of experiment.
'ChIP-seq', 'ChEC-seq', or 'SChEC-seq'.}

\item{paired}{Paired-end status of the reads.
Either TRUE or FALSE.}

\item{ncores}{The number of cores/threads to use.}

\item{genome_dir}{The directory of the genome index.}

\item{genome_index}{The location of the bowtie2 genome index}

\item{genome_assembly}{The directory and file name of the
the gnome assembly FASTA file.}

\item{alignment_dir}{The directory containing the aligned reads.}

\item{peak_dir}{The directory containing the called peaks.}

\item{compare}{For deeptools, do you want to use bamCompare?}

\item{plot_dir}{directory for plotting outputs}

\item{comp_op}{Single operation for comparison of coverage}
}
\description{
Set Settings.
}
