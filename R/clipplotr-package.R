#' clipplotr: comparative visualisation of CLIP data
#'
#' Load CLIP crosslink, auxiliary (e.g. peak) and coverage tracks once with
#' [load_tracks()], then plot any number of regions with [plot_region()].
#' The command-line tool (`./clipplotr`) is a thin wrapper around these
#' functions; see [run_cli()].
#'
#' @keywords internal
#' @import GenomicRanges IRanges S4Vectors GenomeInfoDb rtracklayer GenomicFeatures ggplot2
#' @importFrom data.table data.table as.data.table rbindlist setnames setkey frollmean := .N fread
#' @importFrom patchwork plot_layout
#' @importFrom cowplot theme_minimal_grid theme_minimal_vgrid
#' @importFrom ggthemes scale_colour_tableau scale_fill_tableau
#' @importFrom AnnotationDbi loadDb saveDb
#' @importFrom stats setNames
#' @importFrom utils globalVariables
"_PACKAGE"

# data.table syntax is used inside package functions
.datatable.aware <- TRUE

# Columns referred to with non-standard evaluation in data.table and ggplot2
utils::globalVariables(c(
  "centre", "end", "exp", "gene", "gene_id", "gene_name", "group", "grp",
  "itemRgb", "libSize", "norm", "sample", "score", "smoothed", "start",
  "transcript_id", "type", "width", "x1", "x2", "y1", "y2", "ymax", "ymin"
))
