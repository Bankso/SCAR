
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
#'
#' @export

make_bgs <- function(
  SCAR_obj,
  outdir = getwd(),
  pair_lr = FALSE,
  frag_size = FALSE
 ) {
  ## Input checks.
  paired_status <- as.logical(pull_setting(SCAR_obj, "paired"))
  
  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  ## Make sure we have the correct BAM directory
  set_settings(SCAR_obj, alignment_dir = "./aligned")
  
  ## Prepare command.
  commands <- imap(samples, function(x, y) {
	command <- str_c(
      "bedtools",
      "genomecov", "-ibam",
	  str_c(pull_setting(SCAR_obj, "alignment_dir"), str_c(x, ".bam")),
	  "-bg", 
      sep = " "
    )
	
	if ((paired_status) && (pair_lr)) {
      command <- str_c(command, "-pc", sep = " ")
    }
	
	if ((paired_status) && (frag_size)) {
      command <- str_c(command, "-fs", sep = " ")
    }
	
    return(command)
  })
  
  ## Run commands.
  print_message("Creating bedgraph files from BAMs.")
  walk(commands, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

 ## Add settings to SCAR object.
  print_message("Assigning alignment dir to outdir")
  SCAR_obj <- set_settings(SCAR_obj, alignment_dir = outdir)

  ## Add bam files to sample_sheet.
  print_message("Assigning bedgraphs to sample sheet")
  SCAR_obj <- add_bgs(SCAR_obj, alignment_dir = outdir)

  ## Return the SCAR object.
  return(SCAR_obj)

}