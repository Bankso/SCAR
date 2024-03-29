
#' Plot matrix output as heatmap
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

plot_heatmap <- function(
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
	  											 	"--plotHeight", 10, "--plotWidth", 10,
	  											 "-T", str_c(y, "signal in given BED regions", sep = " "),
	  											 "--colorMap", "cividis",
	  											 "--boxAroundHeatmaps", "no",
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
  											 	str_c("--hclust", clust)
  											 },
  											 "--plotHeight", 10, "--plotWidth", 10,
  											 "-T", "Sample signal in given BED regions",
  											 "--colorMap", "cividis",
  											 "--boxAroundHeatmaps", "no",
  											 sep = " ")
  		}
  	}
	print_message("Deeptools - building heatmap from matrix")
	system2("plotHeatmap", args=command, stderr=str_c(outdir, y, "_log.txt"))
  	}
  )

  ## Return SCAR object.
  return(SCAR_obj)

}
