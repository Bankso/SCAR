
#' Bowtie2 Index Generation
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory for genome index.
#' @param genome_assembly Path to genome fasta file.
#' @param index_name The naming structure given to the index.
#'
#' @export

bowtie2_index <- function(
  SCAR_obj,
  genome_assembly,
  outdir = getwd(),
  index_name = "bowtie2_index"
) {

  ## Make sure output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare Bowtie2 index command.
  command <- str_c(
    "bowtie2-build",
    "-f", genome_assembly,
    "--threads", pull_setting(SCAR_obj, "ncores"),
    str_c(outdir, index_name),
    sep = " "
  )

  ## Run the command.
  print_message("Creating the Bowtie2 genome index.")
  system(command)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Store the genome directory.
  SCAR_obj <- set_settings(
  	SCAR_obj,
  	genome_dir = outdir,
  	genome_index = str_c(outdir, index_name),
  	genome_assembly = genome_assembly
  	)

  print_message("Index location stored in settings.")

  ## Return the SCAR object.
  return(SCAR_obj)

}

#' Bowtie2 Alignment
#'
#' @importFrom purrr walk imap map
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory for aligned reads.
#' @param alignment_mode Either 'end-to-end' or 'local'.
#' @param min_fragment Minimum fragment length (paired end).
#' @param max_fragment Maximum fragment length (paired end).
#' @param mapq_val Minimum MAPQ value to accept reads. Used to filter via
#'                 samtools post alignment
#' @param max_memory Maximum memory per thread for samtools.
#'
#' @export

bowtie2_align <- function(
  SCAR_obj,
  outdir = getwd(),
  alignment_mode = "end-to-end",
  min_fragment = NA,
  max_fragment = NA,
  mapq_val = 10,
  max_memory = "1G"
) {

  ## Input checks.
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))

  ## Create output directory if it exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  if (paired_status) {
    samples <- split(
      SCAR_obj@sample_sheet[, .(sample_name, file_1, file_2)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    if (any(!is.na(SCAR_obj@sample_sheet[["control_file_1"]]))) {
      controls <- split(
        unique(SCAR_obj@sample_sheet[
          !is.na(control_file_1),
          .(control_name, control_file_1, control_file_2)
        ]),
        by = "control_name",
        keep.by = FALSE
      )
      controls <- map(controls, as.character)
      samples <- c(samples, controls)
    }
  }

  else {
    samples <- split(
      SCAR_obj@sample_sheet[, .(sample_name, file_1)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    if (any(!is.na(SCAR_obj@sample_sheet[["control_file_1"]]))) {
      controls <- split(
        unique(SCAR_obj@sample_sheet[
          !is.na(control_file_1),
          .(control_name, control_file_1)
        ]),
        by = "control_name",
        keep.by = FALSE
      )
      controls <- map(controls, as.character)
      samples <- c(samples, controls)
    }
  }

  ## Prepare bowtie2 alignment command.

  print_message("Aligning the FASTQ reads to the genome using Bowtie2")
  iwalk(samples, function(x, y) {
    command <- str_c(
      "-x", pull_setting(SCAR_obj, "genome_index"),
      "-S", str_c(outdir, y, ".sam"),
      "--phred33",
      "--no-unal",
      "-p", pull_setting(SCAR_obj, "ncores"),
      sep = " "
    )

    if (paired_status) {
      command <- str_c(
        command,
        "--no-mixed",
        "--no-discordant",
        "-1", x[1],
        "-2", x[2],
        sep = " "
      )

      if (!is.na(min_fragment)) {
        command <- str_c(command, "-I", min_fragment, sep = " ")
      }
      if (!is.na(max_fragment)) {
        command <- str_c(command, "-X", max_fragment, sep = " ")
      }
    } else {
      command <- str_c(command, "-U", x, sep = " ")
    }

    system2(
      "bowtie2",
      args = command,
      stderr = str_c(outdir, y, "_log.txt")
    )
  })

  ## Make coordinate sorted and indexed bams.
  print_message("Coordinate sorting, filtering and indexing the BAMs")
  walk(names(samples), function(x) {
    command <- str_c(
      "samtools", "sort",
      "-m", max_memory,
      "-@", pull_setting(SCAR_obj, "ncores"),
      "-o", str_c(outdir, str_c(x, ".bam")),
      "-O", "BAM",
      str_c(outdir, str_c(x, ".sam")),
      sep = " "
    	)
    	system(command)

    command <- str_c(
      "samtools", "view", "-q", mapq_val,
      "-b",
      str_c(outdir, str_c(x, ".bam")),
      ">", str_c(outdir, str_c(x, "_sorted.bam")),
      sep = " "
  		)
    	system(command)

    command <- str_c(
      "samtools", "index",
      str_c(outdir, str_c(x, "_sorted.bam")),
      sep = " "
    	)
    	system(command)
  })

  ## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning bams to sample sheet")
  SCAR_obj <- add_bams(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)

}
