
#' Compare two bigWig input files
#' @importFrom purrr map walk pwalk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.'
#' @param comp_op Operation for bamCompare
#' @param bin_size Bin size for coverage summary.
#' @param skip_zeros Skip bins where both files have 0 coverage, TRUE/FALSE
#' @param skip_non_covered Should regions without coverage be skipped.
#' @param scale_factors Takes a named vector, with the name being the
#'   sample name, and the value being the scale factor.
#' @param out_type Output file format, bedgraph or bigwig
#' @param roi A specific region of interest to limit the process to
#'
#' @export

compare_bws <- function(
  SCAR_obj,
  outdir = getwd(),
  comp_op = NA,
  bin_size = 1,
  skip_zeros = TRUE,
  skip_non_covered = FALSE,
  scale_factors = NA,
  out_type = NA,
  roi = NA
) {

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams
  samples <- split(
			SCAR_obj@sample_sheet[, .(sample_name, sample_cov, control_cov)],
			by = "sample_name",
			keep.by = FALSE
			)

  samples <- map(samples, as.character)


  ## Prepare command.
  iwalk(samples, function(x, y) {
	    command <- str_c(
					"-b1", x[1],
					"-b2", x[2],
					"--operation", comp_op,
					"-o", str_c(
	    			outdir, y, "_", comp_op, "_control.cov", sep = ""),
					"-bs", bin_size,
					"-of", out_type,
					"-p", pull_setting(SCAR_obj, "ncores"),
					sep = " ")

  	if (all(!is.na(scale_factors))) {
      command <- str_c(command, str_c(
        "--scaleFactors", scale_factors[y], sep = " "), sep = " ")
			}

  	if (skip_zeros == TRUE) {
      command <- str_c(
        command, "--skipZeroOverZero", sep = " ")
			}

    if (skip_non_covered == TRUE) {
      command <- str_c(
        command, "--skipNonCoveredRegions", sep = " ")
			}

	if (!is.na(roi)) {
      command <- str_c(
        command, "-r", roi, sep = " ")
			}
	
	print_message(command)
	print_message("Deeptools - building comparison tracks from input bigWigs")
	system2("bigwigCompare", args=command, stderr=str_c(outdir, y, "_log.txt"))

  }
  )

  ## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  SCAR_obj <- set_settings(SCAR_obj, compare = TRUE)

  ## Add new comparison tracks to sample_sheet.
  print_message("Assigning comparison tracks to sample sheet")
  SCAR_obj <- add_cov(SCAR_obj, alignment_dir = outdir)

  ## Return SCAR object.
  return(SCAR_obj)

}
