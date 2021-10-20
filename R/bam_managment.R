#' Add BAMs to Sample Sheet
#' @import stringr
#'
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory of aligned reads.
#'
#' @export

add_bams <- function(
  SCAR_obj,
  alignment_dir
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  
  if (!str_detect(alignment_dir, "/$")) {
  alignment_dir <- str_c(alignment_dir, "/")
  }
  
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  ## Add BAM
  
  if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq")) {
    sample_sheet[
		, sample_bams := str_c(alignment_dir, sample_name, ".bam")
	]
    
	sample_sheet[
		!is.na(control_file_1), control_bams := str_c(alignment_dir, control_name, ".bam")
	]
  }

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}