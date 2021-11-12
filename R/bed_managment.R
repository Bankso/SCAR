
#' Add BEDs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param peak_dir Peak directory for SEACR output
#'
#' @export

add_beds <- function(
  SCAR_obj,
  peak_dir
) {

  ## Grab some info from object and prepare inputs.
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  sample_sheet[,
    sample_beds := str_c(peak_dir, sample_name, ".bed")]

  ## Add new sample sheet back to zent object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(SCAR_obj)
}
