
#' Make bedgraphs
#'
#' @importFrom purrr imap
#' @importFrom purrr walk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param pair_lr TRUE or FALSE, calculate coverage as number of
#'        frags covering each bp (paired)
#' @param frag_size Use stated fragment size from pairs instead
#'        of read length (paired)
#' @param chrom_file bedtools genome file path
#'
#' @export

make_bgs <- function(
  SCAR_obj,
  outdir = getwd(),
  pair_lr = FALSE,
  frag_size = FALSE,
  chrom_file = chrom_file
	)
{
  ## Input checks.
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams.
  samples <- split(
  	SCAR_obj@sample_sheet[, .(sample_name, sample_bams)],
  	by = "sample_name",
  	keep.by = FALSE
  	)
  samples <- map(samples, as.character)

  controls <- split(
  	unique(SCAR_obj@sample_sheet[, .(control_name, control_bams)]
  				 ),
  		by = "control_name",
  		keep.by = FALSE
  	)
  	controls <- map(controls, as.character)
  	samples <- c(samples, controls)

  ## Conversion of bams to bedgraphs
  iwalk(samples, function(x, y) {

  	## Step 1: convert bam to bed
  	command <- str_c(
  		"bamtobed",
  		"-bedpe",
  		"-i",
  		x,
  		">", str_c(outdir, str_c(y, ".bed")),
  		sep = " "
  		)
  	system2("bedtools", args=command, stderr=str_c(outdir, y, "_log.txt"))


  	## Step 2: remove extra info from bed file
  	command <- str_c(
  		"'$1==$4 && $6-$2 < 1000 {print $0}'",
  		str_c(outdir, str_c(y, ".bed")),
  		">", str_c(outdir, str_c(y, ".clean.bed")),
  		sep = " "
  		)
  	system2("awk", args=command, stderr=str_c(outdir, y, "_log.txt"))


  	## Step 3: convert to fragments
  	command <- str_c(
  		"cut",
  		"-f",
  		"1,2,6",
  		str_c(outdir, str_c(y, ".clean.bed")),
  		"|",
  		"sort",
  		"-k1,1",
  		"-k2,2n",
  		"-k3,3n",
  		">", str_c(outdir, str_c(y, ".fragments.bed")),
  		sep = " "
  		)
  	system(command)


  	## Step 4: convert to bedgraph via bedtools
  	command <- str_c(
  		"genomecov",
			"-bg",
			"-i",
			str_c(outdir, str_c(y, ".fragments.bed")),
			"-g",
			chrom_file,
			">", str_c(outdir, str_c(y, ".fragments.bedgraph")),
    	sep = " "
    	)

		print_message("bedtools - converting bams to bedgraphs for SEACR")
		system2("bedtools", args=command, stderr=str_c(outdir, y, "_log.txt"))
  })

	## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning bedgraphs to sample sheet")
  SCAR_obj <- add_bgs(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)
  }
