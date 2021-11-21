#' Protected fragment analysis
#'
#' @importFrom purrr imap
#' @importFrom purrr walk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param in_bam path to an input bam, leave blank if none
#' @param in_bg path to an input bedgraph, leave blank if none
#' @param bed_file A bed file of regions to be checked for local low signal sites in aligned bams
#'
#' @export

pf_analysis <- function(
  SCAR_obj,
  outdir = getwd(),
  in_bam = NA,
  in_bg = NA,
  bed_file = bed_file
	)
{

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	if (!is.na(in_bam) && !is.na(in_bg)) {
		print_message("bam and bg detected as inputs - not allowed - exiting")
		command <- "exit"
		system(command)
	}

  ## Conversion of bams to bedgraphs
  if (!is.na(in_bam)) {
		command <- str_c(
			"genomecov",
			"-ibam",
			in_bam,
			"-bga|awk",
			"'$4==0'",
			">", str_c(outdir, "in.bg"),
			sep = " "
			)
	print_message("bedtools - processing input bam")
	system2("bedtools", args=command, stderr=str_c(outdir, "log.txt"))
	}

	if (!is.na(in_bg)) {
		command <- str_c(
			"awk",
			"'$4==0'",
			in_bg,
			">", str_c(outdir, "in.bg"),
			sep = " "
			)
	print_message("Processing input bedgraph")
	system(command)
	}


  ## Compare bg to input bed file
  command <- str_c(
  	"intersect",
  	"-a",
  	str_c(outdir, "in.bg"),
  	"-b",
		bed_file,
  	">", str_c(outdir, "bg_lows_in_peaks.bed"),
  	sep = " "
  	)

	print_message("bedtools - finding low bam signal within bed file regions")
	system2("bedtools", args=command, stderr=str_c(outdir, "intersect_log.txt"))

	## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning overlap beds to sample sheet")
  SCAR_obj <- add_ovbed(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)
  }
