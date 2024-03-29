
#' Set Settings.
#'
#' @param SCAR_obj SCAR object, holds data for use in all the functions.
#' @param analysis_type Type of experiment.
#' 'ChIP-seq', 'ChEC-seq', or 'SChEC-seq'.
#' @param paired Paired-end status of the reads.
#'   Either TRUE or FALSE.
#' @param ncores The number of cores/threads to use.
#' @param genome_dir The directory of the genome index.
#' @param genome_index The location of the bowtie2 genome index
#' @param genome_assembly The directory and file name of the
#'   the gnome assembly FASTA file.
#' @param alignment_dir The directory containing the aligned reads.
#' @param peak_dir The directory containing the called peaks.
#' @param compare For deeptools, do you want to use bamCompare?
#' @param plot_dir directory for plotting outputs
#' @param comp_op Single operation for comparison of coverage
#'
#' @export

set_settings <- function(
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
) {

  settings <- copy(SCAR_obj@settings)

  if (!is.na(analysis_type)) {
    settings[parameter == "analysis_type", value := analysis_type]
  }

  if (!is.na(paired)) {
    settings[parameter == "paired", value := paired]
  }

  if (!is.na(ncores)) {
    settings[parameter == "ncores", value := ncores]
  }

  if (!is.na(genome_dir)) {
    settings[parameter == "genome_dir", value := genome_dir]
  }

  if (!is.na(genome_index)) {
    settings[parameter == "genome_index", value := genome_index]
  }

  if (!is.na(genome_assembly)) {
    settings[parameter == "genome_assembly", value := genome_assembly]
  }

  if (!is.na(alignment_dir)) {
    settings[parameter == "alignment_dir", value := alignment_dir]
  }

  if (!is.na(peak_dir)) {
    settings[parameter == "peak_dir", value := peak_dir]
  }

  if (!is.na(compare)) {
    settings[parameter == "compare", value := compare]
  }

  if (!is.na(plot_dir)) {
  	settings[parameter == "plot_dir", value := plot_dir]
  }
  
  if (!is.na(comp_op)) {
  	settings[parameter == "comp_op", value := comp_op]
  }

  SCAR_obj@settings <- settings
  return(SCAR_obj)
}

#' Pull Settings
#'
#' @param SCAR_obj SCAR object, holds data for use in all the functions.
#' @param setting Setting to pull.
#'
#' @export

pull_setting <- function(
  SCAR_obj,
  setting
) {

  setting_value <- SCAR_obj@settings[parameter == setting, value]
  return(setting_value)

}
