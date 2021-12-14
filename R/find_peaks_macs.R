
#' MACS peak calling
#' @importFrom purrr map walk pwalk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param process Name of the MACS function you want, def. is 'callpeaks'
#' @param in_form Format of input, default is to auto-determine
#' @param genome_size Size of mappable genome
#' @param tag_size size of adaptors, default is to auto-determine
#' @param width Width to use for model building scan, optional
#' @param mfold upper and lower limit for model building, optional
#' @param q_val min q value/FDR cutoff, optional
#' @param p_val Set p val instead of q, optional
#' @param broad Peaks are expected to be broad, TRUE or FALSE, default=FALSE
#'
#' @export

call_peaks_macs <- function(
	SCAR_obj,
	outdir = getwd(),
	process = "callpeak",
	in_form = NA,
	genome_size = NA,
	tag_size = NA,
	width = NA,
	mfold = NA,
	q_val = NA,
	p_val = NA,
	broad = FALSE
) {
	
	peak_dir <- str_c(outdir, "macs/")
	
	## Input checks.
	paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))
	analysis_type <- pull_setting(SCAR_obj, "analysis_type")
	
	## Make output directory if it doesn't exist.
	if (!dir.exists(peak_dir)) dir.create(peak_dir, recursive = TRUE)
	
	## Get bams
	samples <- split(
		SCAR_obj@sample_sheet[, .(sample_name, sample_bams, control_bams)],
		by = "sample_name",
		keep.by = FALSE
	)
	samples <- map(samples, as.character)
	
	## Prepare command.
	iwalk(samples, function(x, y) {
		command <- str_c(
			process,
			"-t", x[1],
			"-c", x[2],
			"--outdir", peak_dir,
			"-n", y,
			"-g", genome_size,
			sep = " ")
		
		if (!is.na(in_form)) {
			command <- str_c(
				command, "-f", in_form, sep = " ")
		}
		
		if (!is.na(tag_size)) {
			command <- str_c(
				command, "-s", tag_size, sep = " ")
		}
		
		if (!is.na(width)) {
			command <- str_c(
				command, "--bw", width, sep = " ")
		}
		
		if (!is.na(mfold)) {
			command <- str_c(
				command, "--mfold", mfold, sep = " ")
		}
		
		if (!is.na(q_val)) {
			command <- str_c(
				command, "-q", q_val, sep = " ")
		}
		
		if (!is.na(p_val)) {
			command <- str_c(
				command, "-p", p_val, sep = " ")
		}
		
		if (broad == TRUE) {
			command <- str_c(
				command, "--broad", sep = " ")
		}
		
		print_message(command)
		print_message("MACS2 - finding peaks in sample alignments")
		
		system2("macs2", args=command, stderr=str_c(outdir, y, "_log.txt"))
		
	}
	)
	
	## Add settings to SCAR object.
	print_message("Assigning peak_dir to macs peak_dir")
	SCAR_obj <- set_settings(SCAR_obj, peak_dir = peak_dir)
	
	## Add bw files to sample_sheet.
	print_message("Assigning macs peaks to sample sheet")
	SCAR_obj <- add_beds(SCAR_obj, peak_dir = peak_dir,
											 peak_type = "macs", stringent = FALSE)
	
	## Return SCAR object.
	return(SCAR_obj)
	
}
