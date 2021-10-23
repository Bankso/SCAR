
#' Add bigwigs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory of aligned reads.
#' @param comp_op Operation for bamCompare
#' @export

add_bws <- function(
  SCAR_obj,
  alignment_dir,
  comp_op = NA
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  if (!str_detect(alignment_dir, "/$")) {
    alignment_dir <- str_c(alignment_dir, "/")
  }
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq") && (compare)) {
    sample_sheet[,
      sample_bws := str_c(alignment_dir, sample_name, 
                    (if (!is.na(comp_op)) {str_c("_", comp_op)}), ".bw")
    ]
  }
  
  else if (!as.logical(compare)) {
    sample_sheet[,
      sample_bws := str_c(alignment_dir, sample_name, 
                    (if (!is.na(comp_op)) {str_c("_", comp_op)}), ".bw")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control_bws := str_c(alignment_dir, control_name, 
                    (if (!is.na(comp_op)) {str_c("_", comp_op)}), ".bw")
    ]
  }

  ## Add new sample sheet back to zent object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(SCAR_obj)
}
