#' clipplotr: comparative visualisation of CLIP data
#'
#' Load CLIP crosslink, auxiliary (e.g. peak) and coverage tracks once with
#' [load_tracks()], then plot any number of regions with [plot_region()].
#' The command-line tool (`./clipplotr`) is a thin wrapper around these
#' functions; see [run_cli()].
#'
#' @keywords internal
#' @import GenomicRanges IRanges S4Vectors GenomeInfoDb rtracklayer GenomicFeatures ggplot2
#' @importFrom data.table data.table as.data.table rbindlist setnames setkey frollmean := .N .I fread
#' @importFrom patchwork plot_layout
#' @importFrom AnnotationDbi loadDb saveDb
#' @importFrom stats setNames
#' @importFrom utils globalVariables
"_PACKAGE"

# data.table syntax is used inside package functions
.datatable.aware <- TRUE

# Columns referred to with non-standard evaluation in data.table and ggplot2
utils::globalVariables(c(
  "canonical", "end", "exp", "gene_id", "group", "grp", "itemRgb", "keep",
  "label", "libSize", "mane", "norm", "position", "sample", "score",
  "smoothed", "start", "transcript_id", "transcript_name", "type", "V1", "width", "y",
  "ymax", "ymin"
))
