
#' Add bigwigs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory of aligned reads.
#' @param comp_op Operation for bamCompare
#' @export

add_bws <- function(
  SCAR_obj,
  alignment_dir=alignment_dir,
  comp_op=comp_op
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  compare <- as.logical(pull_setting(SCAR_obj, "compare"))
  
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  if (compare == TRUE) {
    sample_sheet[,
      sample_bws := str_c(
        alignment_dir, sample_name, "_", comp_op, "_control.bw")
    ]
  }
  
  else {
    sample_sheet[,
      sample_bws := str_c(alignment_dir, sample_name, ".bw")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control_bws := str_c(alignment_dir, control_name, ".bw")
    ]
  }

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}
