
#' Retrieve Reads from SRA
#'
#' @importFrom stringr str_c str_detect str_glue
#' @importFrom purrr walk map
#'
#' @param SCAR_obj SCAR object, holds data for use in all the functions..
#' @param outdir Directory to download files to.
#'
#' @export

retrieve_reads <- function(
  SCAR_obj,
  outdir = getwd()
) {

  ## Input checks.
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Make sure outdir exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get vector of accession numbers.
  accessions <- SCAR_obj@sample_sheet[["file_1"]]

  if (analysis_type %in% c("ChIP-seq", "SChEC-seq", "ChEC-seq")) {
    controls <- unique(SCAR_obj@sample_sheet[["control_file_1"]])
    if (any(!is.na(controls))) {
      controls <- controls[!is.na(controls)]
      accessions <- c(accessions, controls)
    }
  }

  ## Prepare command to get data from ENA.
  # command <- str_c(
  #   "enaDataGet",
  #   "-f", "fastq",
  #   "-d", outdir,
  #   sep = " "
  # )
  
  cores <- pull_setting(SCAR_obj, "ncores")
  command <- str_glue(
    "parallel -I ,, -j {cores} enaDataGet -f fastq -d {outdir} ,, ::: {str_c(accessions, collapse=' ')}"
  )

  ## Run the enaDataGet command.
  print_message("Retrieving the FASTQ files.")
  system(command)

  ## Update the sample sheet.
  sample_sheet <- copy(zent_obj@sample_sheet)

  if (paired_status) {
    sample_sheet[, c("file_1", "file_2") := list(
      str_c(outdir, file_1, "/", file_1, "_1.fastq"),
      str_c(outdir, file_1, "/", file_1, "_2.fastq")
    )]

    if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
      sample_sheet[
        !is.na(control_file_1),
        c("control_file_1", "control_file_2") := list(
          str_c(outdir, control_file_1, "/", control_file_1, "_1.fastq"),
          str_c(outdir, control_file_1, "/", control_file_1, "_2.fastq")
        )
      ]
    }
  } else {
    sample_sheet[, file_1 := str_c(outdir, file_1, "/", file_1, ".fastq")]

    if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
      sample_sheet[
        !is.na(control_file_1),
        control_file_1 := str_c(outdir, control_file_1, "/", control_file_1, ".fastq")
      ]
    }
  }

  zent_obj@sample_sheet <- sample_sheet

  ## Unzip the sequences.
  sequences <- zent_obj@sample_sheet[["file_1"]]
  
  if (paired_status) {
    sequences <- c(sequences, zent_obj@sample_sheet[["file_2"]])
  }

  if (
    analysis_type %in% c("ChIP-seq", "ChEC-seq") &&
    any(!is.na(zent_obj@sample_sheet[["control_file_1"]]))
  ) {

    controls_1 <- unique(zent_obj@sample_sheet[["control_file_1"]])
    controls_1 <- controls_1[!is.na(controls_1)]
 
    sequences <- c(sequences, controls_1)

    if (paired_status) {
      controls_2 <- unique(unique(zent_obj@sample_sheet[["control_file_2"]]))
      controls_2 <- controls_2[!is.na(controls_2)]

      sequences <- c(sequences, controls_2)
    }
  }

  print_message("Unzipping the FASTQ files.")
  walk(sequences, function(x) {
    command <- str_c("gunzip", str_c(x, ".gz"), sep = " ")
    system(command)#, ignore.stdout = TRUE, ignore.stderr = TRUE) 
  })

  ## Return the zent object.
  return(zent_obj)

}
