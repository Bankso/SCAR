#' Protected fragment analysis
#'
#' @importFrom purrr imap
#' @importFrom purrr walk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param bed_file A bed file of regions to be checked for local low signal sites in aligned bams
#'
#' @export

pf_analysis <- function(
  SCAR_obj,
  outdir = getwd(),
  bed_file = bed_file
	)
{

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  samples <- split(
      SCAR_obj@sample_sheet[, .(sample_name, sample_bams)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    if (any(!is.na(SCAR_obj@sample_sheet[["control_bams"]]))) {
      controls <- split(
        unique(SCAR_obj@sample_sheet[
          !is.na(control_file_1),
          .(control_name, control_bams)
        ]),
        by = "control_name",
        keep.by = FALSE
      )
      controls <- map(controls, as.character)
      samples <- c(samples, controls)
	}
  ## Conversion of bams to bedgraphs
  iwalk(samples, function(x, y) {

  	command <- str_c(
  		"genomecov",
		"-ibam",
		x,
		"-bga|awk",
		"'$4==0'",
		">", str_c(outdir, str_c(y, "_bam.bg")),
    	sep = " "
    	)

		print_message("bedtools - converting bam to bedgraph for protected sequence analysis")
		system2("bedtools", args=command, stderr=str_c(outdir, y, "_log.txt"))

  	## Compare bg to input bed file
  	command <- str_c(
  		"intersect",
  		"-a",
  		str_c(outdir, str_c(y, "_bam.bg")),
  		"-b",
		bed_file,
  		">", str_c(outdir, str_c(y, "_overlap.bed")),
  		sep = " "
  		)

		print_message("bedtools - finding low bam signal within bed file regions")
		system2("bedtools", args=command, stderr=str_c(outdir, y, "_log.txt"))

  })
	## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning overlap beds to sample sheet")
  SCAR_obj <- add_ovbed(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)
  }
