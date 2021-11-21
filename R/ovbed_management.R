#' Add overlap beds to Sample Sheet
#' @import stringr
#'
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory for files to be dumped to.
#'
#' @export

add_ovbed <- function(
  SCAR_obj,
  alignment_dir
) {

  ## Grab some info from object and prepare inputs.

  sample_sheet <- copy(SCAR_obj@sample_sheet)

  sample_sheet[,
    	ov_bed := str_c(alignment_dir, "bg_lows_in_peaks.bed")]

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}