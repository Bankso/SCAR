#' Add Bedgraphs to Sample Sheet
#' @import stringr
#' 
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory for files to be dumped to.
#'
#' @export

add_bgs <- function(
  SCAR_obj,
  alignment_dir
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  if (!str_detect(alignment_dir, "/$")) {
    alignment_dir <- str_c(alignment_dir, "/")
  }
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  ## Add if RNA-seq experiment.
  if (analysis_type == "RNA-seq") {
    sample_sheet[, bg_files := str_c(
      alignment_dir, sample_name, "_sorted.bedgraph"
    )]
  } else if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq")) {
    sample_sheet[, sample_bgs := str_c(alignment_dir, sample_name, ".bedgraph")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control__bgs := str_c(alignment_dir, control_name, ".bedgraph")
    ]
  }

  ## Add new sample sheet back to zent object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(SCAR_obj)
}

