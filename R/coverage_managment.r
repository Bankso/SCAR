
#' Add coverage tracks to Sample Sheet
#' @import stringr
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory of aligned reads.
#' @param comp_op Operation for bamCompare
#' @export

add_cov <- function(
  SCAR_obj,
  alignment_dir,
  comp_op = 'comp'
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  compare <- as.logical(pull_setting(SCAR_obj, "compare"))
  alignment_dir <- pull_setting(SCAR_obj, 'alignment_dir')
  
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  if (compare == TRUE) {
    sample_sheet[, compared_cov := str_c(
      alignment_dir, sample_name, "_", comp_op, "_control.cov")
      ]
    }
  
  else {
    sample_sheet[,
      sample_cov := str_c(alignment_dir, sample_name, ".cov")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control_cov := str_c(alignment_dir, control_name, ".cov")
    ]
  }

  ## Add new sample sheet back to SCAR object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the SCAR object.
  return(SCAR_obj)
}
