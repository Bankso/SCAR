#' Annotate Peaks
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom rtracklayer import
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom purrr iwalk
#'
#' @param SCAR_obj SCAR object.
#' @param outdir Directory to output the annotated peaks.
#' @param genome_annotation Genome GTF annotation file.
#' @param promoter_downstream Bases downstream of TSS to define promoter.
#' @param promoter_upstream Bases upstream of TSS to define promoter.
#' @param feature_type Either 'gene' or 'transcript'.
#'
#' @export

annotate_peaks <- function(
  SCAR_obj,
  outdir = getwd(),
  genome_annotation,
  promoter_downstream = 1000,
  promoter_upstream = 1000,
  feature_type = "gene"
) {

  ## Make sure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) 

  ## Import the genome annotation file.
  genome_txdb <- makeTxDbFromGFF(genome_annotation)

  ## Import the peak files.
  peak_files <- str_c(
    pull_setting(SCAR_obj, "peak_dir"),
    SCAR_obj@sample_sheet[["sample_name"]],
    "_peaks.narrowPeak"
  )
  names(peak_files) <- SCAR_obj@sample_sheet[["sample_name"]]

  narrowpeak_cols  <-  c(
    signal.value = "numeric",
    p.value.negLog10 = "numeric",
    q.value.negLog10 = "numeric",
    peak = "integer"
  )

  print_message("Annotating the called peaks.")
  iwalk(peak_files, function(x, y) {
    peaks <- import(x, format = "BED", extraCols = narrowpeak_cols)

    annotated_peaks <- annotatePeak(
      peaks,
      tssRegion = c(-promoter_upstream, promoter_downstream),
      TxDb = genome_txdb,
      level = feature_type
    )
    annotated_peaks <- as.data.table(annotated_peaks)

    fwrite(
      annotated_peaks, str_c(outdir, y, "_peaks_annotated.tsv"),
      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
    )
  })

  ## Return the SCAR object.
  return(SCAR_obj)

}
