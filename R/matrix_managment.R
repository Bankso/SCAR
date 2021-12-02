
#' Add matrix to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param plot_dir directory for matrix output
#'
#' @export

add_matrix <- function(
  SCAR_obj,
  plot_dir=plot_dir
) {

	sample_sheet <- copy(SCAR_obj@sample_sheet)

	plot_dir <- pull_setting(SCAR_obj, "plot_dir")

  sample_sheet[,
    sample_matrix := str_c(plot_dir, sample_name, ".matrix")]

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}
