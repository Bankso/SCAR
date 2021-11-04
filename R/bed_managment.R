
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
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq")) {
    sample_sheet[,
      sample_beds := str_c(peak_dir, sample_name, ".bed")
    ]
  }

  ## Add new sample sheet back to zent object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(SCAR_obj)
}
