
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
  sep = ""
  )
{

  ## Make sure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

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
			num_thresh,
			sep = " "
    )

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
	}
  )

  ## Run the commands.
  print_message("Calling peaks from the aligned reads.")
  walk(commands, system2(
  	"bash", args=commands, stderr=str_c(outdir, "SEACR_log.txt"))
  	)

  ## Save the peak directory.
  SCAR_obj <- set_settings(SCAR_obj, peak_dir=outdir)

  ## Add new BEDs to peak_dir
  SCAR_obj <- add_beds(SCAR_obj, peak_dir=outdir, stringent=stringent)

  ## Return the SCAR object.
  return(SCAR_obj)
}
