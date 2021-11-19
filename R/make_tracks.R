
#' Generate coverage tracks
#' @importFrom purrr map walk pwalk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param comp_op Operation for bamCompare
#' @param bin_size Bin size for coverage summary.
#' @param normalize_using Either 'CPM' or 'RPGC'/'RPKM'.
#' @param genome_size Effective genome size, req'd for RPGC/RPKM.
#' @param skip_non_covered Should regions without coverage be skipped.
#' @param min_fragment Minimum fragment length.
#' @param max_fragment Maximum fragment length.
#' @param extend_reads Distance to extend single-end reads.
#'   Set to NA to not extend reads.
#' @param scale_factors Takes a named vector, with the name being the
#'   sample name, and the value being the scale factor.
#'   If set will override 'normalize_using' option.
#' @param split_strands For RNA-seq, whether to split the strands into
#'   positive and minus strand files.
#' @param library_type If split_strands is TRUE, specify library chemistry
#'   as either 'dUTP' or 'ligation'.
#' @param center_reads Should reads be centered based on fragment length? Helps with signal around enriched regions
#' @param out_type Output file format, bedgraph or bigwig
#' @param temp_dir Temporary directory to write files to.
#'
#' @export

make_tracks <- function(
  SCAR_obj,
  outdir = getwd(),
  comp_op = "log2",
  bin_size = 1,
  normalize_using = NA,
  genome_size = NA,
  skip_non_covered = NA,
  min_fragment = NA,
  max_fragment = NA,
  extend_reads = NA,
  scale_factors = NA,
  split_strands = NA,
  library_type = NA,
  center_reads = NA,
  out_type = NA,
  temp_dir = "./temp"
) {

  ## Input checks.
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  compare <- as.logical(pull_setting(SCAR_obj, "compare"))

  ## Set temporary directory.
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  Sys.setenv(TMPDIR=temp_dir)

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams.
  samples <- split(
  	SCAR_obj@sample_sheet[, .(sample_name, sample_bams, control_bams)],
  	by = "sample_name",
  	keep.by = FALSE
  	)
   samples <- map(samples, as.character)


  ## Prepare command.
  if (compare == TRUE){
  iwalk(samples, function(x, y) {
	    command <- str_c(
	    "-b1", x[1],
	    "-b2", x[2],
	    "--operation", comp_op,
	    "-bs", bin_size,
		"-of", out_type,
	    "-o", str_c(
	    	outdir, y, "_", comp_op, "_control.cov", sep = ""),
	          "-p", pull_setting(SCAR_obj, "ncores"), sep = " ")

	if (compare) {
      command <- str_c(
        command, "--scaleFactorsMethod", "None", sep = " ")
    }

	if (all(is.na(scale_factors)) && !is.na(normalize_using)) {
      command <- str_c(
        command, "--normalizeUsing", normalize_using, sep = " ")
    }

    if (all(!is.na(scale_factors))) {
      command <- str_c(command, str_c(
        "--scaleFactor", scale_factors[y], sep = " "), sep = " ")
    }

    if (!is.na(genome_size)) {
      command <- str_c(
        command, "--effectiveGenomeSize", genome_size, sep = " ")
    }

    if (skip_non_covered) {
      command <- str_c(
        command, "--skipNonCoveredRegions", sep = " ")
    }

	if (center_reads) {
      command <- str_c(
        command, "--centerReads", sep = " ")
    }

    if (!is.na(min_fragment)) {
      command <- str_c(
        command, "--minFragmentLength", min_fragment, sep = " ")
    }
    if (!is.na(max_fragment)) {
      command <- str_c(
        command, "--maxFragmentLength", max_fragment, sep = " ")
    }

    if (!is.na(extend_reads) && paired_status) {
      command <- str_c(command, "-e", sep = " ")
    } else if (!is.na(extend_reads) && !paired_status) {
      command <- str_c(command, "-e", extend_reads, sep = " ")
    }

	  print_message("Deeptools - building comparison tracks from aligned reads")
	  system2("bamCompare", args=command, stderr=str_c(outdir, y, "_log.txt"))
	}
  )
}

  ## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bw files to sample_sheet.
  print_message("Assigning bws to sample sheet")
  SCAR_obj <- add_cov(SCAR_obj, alignment_dir = outdir)

  ## Return SCAR object.
  return(SCAR_obj)

}
