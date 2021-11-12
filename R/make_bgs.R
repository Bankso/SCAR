
#' Make bedgraphs
#'
#' @importFrom purrr imap
#' @importFrom purrr walk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param pair_lr TRUE or FALSE, calculate coverage as number of
#'        frags covering each bp (paired)
#' @param frag_size Use stated fragment size from pairs instead
#'        of read length (paired)
#'
#' @export

make_bgs <- function(
  SCAR_obj,
  outdir = getwd(),
  pair_lr = FALSE,
  frag_size = FALSE
	)
{
  ## Input checks.
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams.
  samples <- split(
  	SCAR_obj@sample_sheet[, .(sample_name, sample_bams)],
  	by = "sample_name",
  	keep.by = FALSE
  	)
  samples <- map(samples, as.character)

  controls <- split(
  	unique(SCAR_obj@sample_sheet[, .(control_name, control_bams)]
  				 ),
  		by = "control_name",
  		keep.by = FALSE
  	)
  	controls <- map(controls, as.character)
  	samples <- c(samples, controls)

  ## Prepare command.
  iwalk(samples, function(x, y) {
  	command <- str_c(
  		"genomecov",
			"-ibam",
			x,
			"-bg",
    	sep = " "
    	)

		if ((paired_status) && (pair_lr)) {
      command <- str_c(command, "-pc", sep = " ")
    	}

		if ((paired_status) && (frag_size)) {
      command <- str_c(command, "-fs", sep = " ")
			}
		print(command)
		print_message("bedtools - converting bams to bedgraphs for SEACR")
		system2("bedtools", args=command, stderr=str_c(outdir, y, "_log.txt"))
  	}
	)

	## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning bedgraphs to sample sheet")
  SCAR_obj <- add_bgs(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)
  }
