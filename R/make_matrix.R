
#' Convert coverage file to matrix
#' @importFrom purrr map walk pwalk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.'
#' @param primary Main command; either 'reference-point' or 'scale-regions'
#' @param s_n_c If desired, select a stored bigWig set for processing;
#'        "s" for sample, "n" for control, "c" for compare
#' @param in_str A list of input strings to be converted to a bash command for
#'        computeMatrix
#' @param regions A bed file defining regions to be plotted; if none entered,
#'        the peak file from peak finding will be used
#' @param peak_type Peak finder used to make input peaks, 'macs' or 'seacr'
#' @export

make_matrix <- function(
  SCAR_obj,
  outdir = getwd(),
  primary = "reference-point",
  s_n_c = NA,
  in_str = NA,
  regions = NA,
  peak_type = NA
) {

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams
  samples <- split(
			SCAR_obj@sample_sheet[
				, .(
				sample_name,
				sample_cov,
				control_cov,
				compared_cov,
				macs_peaks,
				seacr_peaks)],
			by = "sample_name",
			keep.by = FALSE
			)

  samples <- map(samples, as.character)

  ## Prepare command.
  iwalk(samples, function(x, y) {
	    command <- str_c(
					
	    		primary, "-S",
					
					if (!is.na(s_n_c)) {
						
						if (s_n_c == "s") {
							x[1]
						}
						
						if (s_n_c == "n") {
							x[2]
						}
						
						if (s_n_c == "c") {
							x[3]
						}
						
						else {
							s_n_c
						}
					},
					"-R", 
					if (!is.na(regions)) {
						regions
						}
						else if (peak_type == "macs") {
							x[4]
							}
						else if (peak_type == "seacr") {
							x[5]
						},
					"-o", str_c(outdir, y, ".matrix"),
					sep = " ")

	  	if (!is.na(in_str)) {
	    	command <- str_c(command, str_c(in_str), sep = " ")
	  	}
	    
	    	else {
	  			command <- str_c(command, "--referencePoint", "center", "-b", 500,
	  								 "-a", 500, "-bs", 1, sep = " ")
	    }
	
	print_message(command)
	print_message("Deeptools - building a plottable matrix from inputs")
	system2("computeMatrix", args=command, stderr=str_c(outdir, y, "_mat_log.txt"))

  }
  )

  ## Add settings to SCAR object.
  print_message("Assigning plot dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, plot_dir = outdir)

  ## Add new comparison tracks to sample_sheet.
  print_message("Assigning matrix to sample sheet")
  SCAR_obj <- add_matrix(SCAR_obj, plot_dir = outdir)

  ## Return SCAR object.
  return(SCAR_obj)

}
