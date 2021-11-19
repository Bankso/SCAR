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
    	sample_ovbeds := str_c(alignment_dir, sample_name, "_overlap.bed")]

  sample_sheet[,
    	control_ovbeds := str_c(alignment_dir, control_name, "_overlap.bed")]

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}