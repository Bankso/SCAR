
#' SCAR Class
#'
#' @importFrom data.table data.table
#' @import stringr
#' @slot sample_sheet Sample sheet.
#' @slot settings Various run settings.
#'
#' @export

setClass(
  "SCAR_obj",
  representation(
    sample_sheet = "data.table",
    settings = "data.table"
  ),
  prototype(
    sample_sheet = data.table(),
    settings = data.table()
  )
)

#' SCAR Object Constructor Function
#'
#' @import methods
#' @import data.table
#'
#' @param sample_sheet Either a data.frame or delimited file.
#' @param analysis_type One of 'ChIP-seq', 'ChEC-seq', or 'SChEC-seq'.
#' @param sep If the sample sheet is a file, this specifies the delimiter.
#' @param paired Whether the run is paired (TRUE) or unpaired (FALSE)
#' @param ncores The number of cores/threads to use.
#' @param compare TRUE/FALSE perform bamCompare instead of Coverage
#' @param comp_op Single operation for comparison of coverage
#'
#' @return A SCAR object.
#'
#' @export

SCAR_maker <- function(
  analysis_type,
  sample_sheet,
  sep = "\t",
  paired = NA,
  ncores = NA,
  compare = NA,
  comp_op = NA
) {

  ## Prepare the sample sheet.
  if (is(sample_sheet, "character")) {
    samples <- fread(sample_sheet, sep = sep)
  } else if (is(sample_sheet, "data.frame")) {
    samples <- as.data.table(sample_sheet)
  }

  samples[is.na(samples)] <- NA_character_

  samples[!is.na(control_name),
      c("control_file_1", "control_file_2", "control_name") := list(
        ifelse(str_detect(control_file_1, "(^$|^NA$)"), NA_character_, control_file_1),
        ifelse(str_detect(control_file_2, "(^$|^NA$)"), NA_character_, control_file_2),
        ifelse(str_detect(control_name, "(^$|^NA$)"), NA_character_, control_name)
      )
    ]


  ## prepare run settings.
  run_settings <- data.table(
     parameter = c(
       "analysis_type", "paired", "ncores", "genome_dir",
       "genome_index", "genome_assembly", "alignment_dir",
       "peak_dir", "compare", "plot_dir", "comp_op"
    ),
     value = c(analysis_type, paired, ncores, rep("", 5), compare, "", "")
  )

  ## Create the SCAR object.
  SCAR_obj <- new(
    "SCAR_obj",
    sample_sheet = samples,
    settings = run_settings
  )

  return(SCAR_obj)

}
