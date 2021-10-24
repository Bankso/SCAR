
#' Peak Calling
#'
#' @importFrom purrr imap
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory for peak files.
#' @param num_thresh P-value thresholding for peaks.
#' @param norm TRUE or FALSE - normalization on/off
#' @param stringent TRUE or FALSE - stringent or relaxed analysis
#' @param sep separator for entries, default here is ""
#'
#' @export

call_peaks_SEACR <- function(
  SCAR_obj,
  outdir = getwd(),
  num_thresh = 0.05,
  norm = TRUE,
  stringent = TRUE,
  sep = ''
  ) 
{

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))

  ## Make sure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get your bgs
  samples <- split(
      SCAR_obj@sample_sheet[, .(sample_name, sample_bgs)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    if (any(!is.na(SCAR_obj@sample_sheet[["control_bgs"]]))) {
      controls <- split(
        unique(SCAR_obj@sample_sheet[
          !is.na(control_bgs),
          .(control_name, control_bgs)
        ]),
        by = "control_name",
        keep.by = FALSE
      )
      controls <- map(controls, as.character)
      
	  samples <- c(samples, controls)
    }
  
  ## Create the peak calling command.
  
  commands <- imap(samples, function(x, y) {
	command <- str_c(
		'bash', 
		'SEACR_1.3.sh',
		str_c(
		  pull_setting(SCAR_obj, "alignment_dir"),
		  str_c(x, ".bedgraph")),
		str_c(
		  pull_setting(SCAR_obj, "alignment_dir"),
		  str_c(y, ".bedgraph")),
		num_thresh,
		sep = " "
    )

    if (norm) {
		command <- str_c(
		command, 'norm', sep = sep
		)
    } 
    
    else { 
		command <- str_c(
		command, 'non', sep = sep 
		)
    }
    
    if (stringent) {
		command <- str_c(
		command, 'stringent', sep = sep
		)
    }
	  else {
		command <- str_c(
		command, 'relaxed', sep = sep 
		)
    }
    
	  command <- str_c(
	  command, str_c(outdir, x), sep = sep 
	  )
    return(command)
	}
  )

  ## Run the commands.
  print_message("Calling peaks from the aligned reads.")
  walk(commands, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Save the peak directory.
  SCAR_obj <- set_settings(SCAR_obj, peak_dir = outdir)
  
  ## Add new BEDs to peak_dir
  SCAR_obj <- add_beds(SCAR_obj, peak_dir = outdir)
  
  ## Return the SCAR object.
  return(SCAR_obj)
}