
#' Add BAMs to Sample Sheet
#'
#' @param SCAR_obj SCAR object.
#' @param alignment_dir Directory of aligned reads.
#'
#' @export

add_bams <- function(
  SCAR_obj,
  alignment_dir=alignment_dir
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  sample_sheet <- copy(SCAR_obj@sample_sheet)

  if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq")) {
    sample_sheet[,
      sample_bams := str_c(alignment_dir, sample_name, "_sorted.bam")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control_bams := str_c(alignment_dir, control_name, "_sorted.bam")
    ]
  }

  ## Add new sample sheet back to zent object.
  SCAR_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(SCAR_obj)
}
