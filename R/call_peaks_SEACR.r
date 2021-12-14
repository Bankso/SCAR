
#' Peak Calling
#'
#' @importFrom purrr imap
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory for peak files.
#' @param num_thresh P-value thresholding for peaks.
#' @param norm TRUE or FALSE - inputs are normalized or not
#' @param stringent TRUE or FALSE - stringent or relaxed analysis
#' @param sep separator for entries, default here is ""
#'
#' @export

call_peaks_SEACR <- function(
  SCAR_obj,
  outdir = getwd(),
  num_thresh = NA,
  norm = NA,
  stringent = NA,
  sep = ""
  )
{
	peak_dir <- str_c(outdir, "seacr/")
	
  ## Make sure the output directory exists.
  if (!dir.exists(peak_dir)) dir.create(peak_dir, recursive = TRUE)

  ## Get your bgs
  samples <- split(
  	SCAR_obj@sample_sheet[, .(sample_name, sample_bgs, control_bgs)],
  	by = "sample_name",
  	keep.by = FALSE
  )
  samples <- map(samples, as.character)


  ## Create the peak calling command.

  commands <- imap(samples, function(x, y) {
	command <- str_c(
		'SEACR_1.3.sh',
		x[1],
		x[2],
		sep = " "
		)

	if (!is.na(num_thresh) && is.na(x[2])) {
		command <- str_c(
			command, num_thresh, sep = " "
			)
		}

    if (norm) {
		command <- str_c(
			command, 'norm', sep = " "
			)
		}

    else {
		command <- str_c(
			command, 'non', sep = " "
			)
		}

    if (stringent) {
		command <- str_c(
			command, 'stringent', sep = " "
			)
		}
	  else {
		command <- str_c(
			command, 'relaxed', sep = " "
			)
		}
	command <- str_c(
			command, str_c(peak_dir, y), sep = " "
			)
	}
  )

  ## Run the commands.
  print_message("SEACR- calling peaks from the aligned reads")
  walk(commands, system2(
  	"bash", args=commands, stderr=str_c(peak_dir, "SEACR_log.txt"))
  	)

  ## Save the peak directory.
  SCAR_obj <- set_settings(SCAR_obj, peak_dir=peak_dir)

  ## Add new BEDs to peak_dir
  SCAR_obj <- add_beds(SCAR_obj, peak_dir = peak_dir,
  										 peak_type = "seacr", stringent = stringent)

  ## Return the SCAR object.
  return(SCAR_obj)
}
