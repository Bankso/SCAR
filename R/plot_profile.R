
#' Plot matrix output in profile
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Output directory.
#' @param matrix If desired, provide a matrix for processing; by default,
#'        the stored matrix is used
#' @param plot_name name for plot with extension (.pdf or .png)
#' @param plot_opts A list of input strings to be converted to a bash command for
#'        the plot function
#' @param clust If given, data will be clustered into the number of groups
#'        equal to the value entered (integers only) via hclust
#'
#' @export

plot_profile <- function(
  SCAR_obj,
  outdir = getwd(),
  matrix = NA,
  plot_name = NA,
  plot_opts = NA,
  clust = NA
) {

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  samples <- split(
			SCAR_obj@sample_sheet[
				, .(
				sample_name,
				sample_matrix)],
			by = "sample_name",
			keep.by = FALSE
  		)

  samples <- map(samples, as.character)

  ## Prepare basic command.
  iwalk(samples, function(x, y) {
  	if (is.na(matrix)) {
	    command <- str_c(
	    	"-m", x,
	    	"-o", str_c(outdir, plot_name),
				"--outFileSortedRegions", str_c(outdir, y, "_sorted.regions"),
					sep = " ")

	  	if (!is.na(plot_opts)) {
	    	command <- str_c(command, str_c(plot_opts), sep = " ")
	    }	else {
	    		## Standard command, used if no other string is provided
	  			command <- str_c(command,
	  											 if (!is.na(clust)) {
	  											 	str_c("--hclust", clust, sep = " ")
	  											 },
	  											 	"--plotHeight", '10', "--plotWidth", '10',
	  											 "-T", 
	  											 str_c(y, "signal in given BED regions", sep = " "),
	  											 sep = " ")
	    										}
  	}
  	if (!is.na(matrix)) {
  		command <- str_c(
  			"-m", matrix,
  			"-o", str_c(outdir, plot_name),
  			"--outFileSortedRegions", str_c(outdir, "plot_sorted.regions"),
  			sep = " ")

  		if (!is.na(plot_opts)) {
  			command <- str_c(command, str_c(plot_opts), sep = " ")
  		}	else {
  			## Standard command, used if no other string is provided
  			command <- str_c(command,
  											 if (!is.na(clust)) {
  											 	str_c("--hclust", clust, sep = " ")
  											 },
  											 "--plotHeight", '10', "--plotWidth", '10',
  											 "-T", "Sample signal in given BED regions",
  											 sep = " ")
  		}
  	}
  
  print_message(command)
	print_message("Deeptools - building a profile plot from inputs")
	system2("plotProfile", args=command, 
					stderr=str_c(outdir, y, "_profile_log.txt"))
  	}
  )

  ## Return SCAR object.
  return(SCAR_obj)

}
