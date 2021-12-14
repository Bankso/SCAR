
#' Add BEDs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param peak_dir Peak directory for SEACR output
#' @param stringent was SEACR run with the 'stringent' flag?
#'
#' @export

add_beds <- function(
  SCAR_obj,
  peak_dir,
  peak_type=peak_type,
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

  if (!is.na(peak_type)) {
  	if (peak_type == "seacr") {
  		sample_sheet[,
    	seacr_peaks := str_c(peak_dir, sample_name, level, ".bed")]
  	}
  	if (peak_type == "macs") {
  		sample_sheet[,
  		macs_peaks := str_c(peak_dir, sample_name, "_summits.bed")]
  	}
  }
  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}
