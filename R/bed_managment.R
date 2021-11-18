
#' Add BEDs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param peak_dir Peak directory for SEACR output
#' @param stringent was SEACR run with the "stringent" flag?
#'
#' @export

add_beds <- function(
  SCAR_obj,
  peak_dir,
  stringent=stringent
) {

  ## Grab some info from object and prepare inputs.
  if (stringent) {
  	level <- ".stringent"
  }
  else {
  	level <- ".relaxed"
  }

	sample_sheet <- copy(SCAR_obj@sample_sheet)

  sample_sheet[,
    sample_beds := str_c(peak_dir, sample_name, level, ".bed")]

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}
