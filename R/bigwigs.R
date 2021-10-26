
#' Generate Bigwigs
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param comp_op Operation for bamCompare
#' @param bin_size Bin size for coverage summary.
#' @param normalize_using Either 'CPM' or 'RPGC'.
#' @param genome_size Effective genome size, req'd for RPGC.
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
#' @param temp_dir Temporary directory to write files to.  
#'
#' @export

make_bigwigs <- function(
  SCAR_obj,
  outdir = getwd(),
  comp_op = "log2",
  bin_size = 1,
  normalize_using = NA,
  genome_size = NA,
  skip_non_covered = TRUE,
  min_fragment = NA,
  max_fragment = NA,
  extend_reads = NA,
  scale_factors = NA,
  split_strands = NA,
  library_type = NA,
  temp_dir = "./temp"
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))
  analysis_type <- pull_setting(SCAR_obj, "analysis_type")
  compare <- as.logical(pull_setting(SCAR_obj, "compare"))

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams.
  if (analysis_type %in% c("ChIP-seq", "ChEC-seq", "SChEC-seq")) {
    samples <- split(
      SCAR_obj@sample_sheet[, .(sample_name, sample_bams)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    if(any(!is.na(SCAR_obj@sample_sheet[["control_bams"]]))) {
      controls <- split(
        unique(SCAR_obj@sample_sheet[
          !is.na(control_bams),
          .(control_name, control_bams)
        ]),
        by = "control_name",
        keep.by = FALSE
      )
      controls <- map(controls, as.character)
      samples <- c(samples, controls)
    }
  } 

  ## Prepare command.
  commands <- iwalk(samples, function(x, y) {
	  command <- (
	    if (compare != TRUE) {
	    print_message("bamCoverage selected based on inputs")
	    str_c("bamCoverage",
	          "-b", x,
	          "-of", "bigwig",
	          "-bs", bin_size,
	          "-o", str_c(outdir, (SCAR_obj@sample_sheet[, .(sample_name)]), 
	                      ".bw", sep = ""),
	          "-p", pull_setting(SCAR_obj, "ncores"),
	          sep = " ")
	  }
	  else {
	    print_message("bamCompare selected based on inputs")
	    str_c("bamCompare",
	          "-b1", x,
	          "-b2", y,
	          "--operation", comp_op,
	          "-bs", bin_size,
	          "-o", str_c(outdir, (SCAR_obj@sample_sheet[
	              , .(sample_name)]), "_", comp_op, "_control.bw", sep = ""),
	          "-p", pull_setting(SCAR_obj, "ncores"),
	          sep = " ")
	  })
	
    
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

    return(command)
  })
  
  ## Set temporary directory.
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  Sys.setenv(TMPDIR=temp_dir)

  ## Run commands.
  print_message("Creating the BIGWIG coverage tracks.")
  walk(commands, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)
  
  ## Add bw files to sample_sheet.
  print_message("Assigning bws to sample sheet")
  SCAR_obj <- add_bws(SCAR_obj, alignment_dir = outdir)
  
  ## Return SCAR object.
  return(SCAR_obj)

}
