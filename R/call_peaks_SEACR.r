
#' Peak Calling
#'
#' @importFrom purrr pmap
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory for peak files.
#' @param num_thresh P-value thresholding for peaks.
#' @param norm TRUE or FALSE - normalization on/off
#' @param stringent TRUE or FALSE - stringent or relaxed analysis
#' @param sep separator for entries, default here is ""
#'
#' @export

call_peaks_SEACR <- function
	(
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

  ## Create the peak calling command.
  commands <- pmap(SCAR_obj@sample_sheet, function(...) {
    args <- list(...)

    command <- str_c(
      'bash', 
      'SEACR_1.3.sh',
	    args$sample_bgs,
      args$control_bgs,
	    num_thresh,
      sep = sep
    )

    if (norm) {
		  command <- str_c(
		  command, 
		  'norm',
		  sep = ,
      )
    } 
    
    else { 
		command <- str_c(command, 'non', sep = sep )
    }
    
    if (stringent) {
      command <- str_c(command, 'stringent', sep = sep )
    } else {
      command <- str_c(command, 'relaxed', sep = sep )
    }

    return(command)
  }
  )

  ## Run the commands.
  print_message("Calling peaks from the aligned reads.")
  walk(commands, system, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Save the peak directory.
  SCAR_obj <- set_settings(SCAR_obj, peak_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)
}