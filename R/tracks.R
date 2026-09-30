#' Load CLIP, auxiliary and coverage tracks
#'
#' Reads all input files once so that many regions can be plotted with
#' [plot_region()] without re-reading them. Crosslink files are read in full
#' (library sizes need the whole file); coverage bigWigs are only read for each
#' region when plotting.
#'
#' @param xlinks Crosslink files: all iCount bedgraphs (`.bedgraph`, score sign
#'   gives the strand) or all BED files, optionally gzip compressed.
#' @param gtf Reference GTF file, or an annotation from [load_annotation()].
#' @param labels,colours,groups Optional label, colour and group for each
#'   crosslink file. Labels default to the file name without extension.
#' @param size_factors Optional size factor for each crosslink file, used by
#'   `custom` and `custom_maxpeak` normalisation.
#' @param auxiliary,auxiliary_labels Optional BED files of features (e.g. peaks)
#'   and their labels.
#' @param coverage,coverage_labels,coverage_colours,coverage_groups Optional
#'   bigWig coverage files (e.g. RNA-seq, on the same strand as the region) and
#'   their labels, colours and groups.
#' @param cache_dir Directory for the cached annotation, see [load_annotation()].
#' @param no_ucsc If `FALSE` (default), chromosome names are converted to UCSC
#'   style. Set to `TRUE` for non-standard sequences, and make sure all files
#'   use the same names.
#' @return A `clipplotr_tracks` object.
#' @export
load_tracks <- function(xlinks, gtf,
                        labels = NULL, colours = NULL, groups = NULL, size_factors = NULL,
                        auxiliary = NULL, auxiliary_labels = NULL,
                        coverage = NULL, coverage_labels = NULL, coverage_colours = NULL, coverage_groups = NULL,
                        cache_dir = NULL, no_ucsc = FALSE) {

  # Check options before reading anything
  if(length(xlinks) == 0) clipplotr_error("No iCLIP bedgraphs supplied")
  check_length(labels, length(xlinks), "labels", "xlinks")
  check_length(colours, length(xlinks), "colours", "xlinks")
  check_length(groups, length(xlinks), "groups", "xlinks")
  check_length(size_factors, length(xlinks), "size_factors", "xlinks")
  if(!is.null(size_factors)) {
    size_factors <- suppressWarnings(as.numeric(size_factors))
    if(any(is.na(size_factors)) || any(size_factors <= 0)) clipplotr_error("size_factors need to be positive numbers")
  }
  if(is.null(labels)) labels <- make.unique(label_from_file(xlinks, "bedgraph|bed"), sep = " ")
  check_unique(labels, "labels")

  check_length(auxiliary_labels, length(auxiliary), "auxiliary_labels", "auxiliary")
  if(length(auxiliary) > 0 && is.null(auxiliary_labels)) auxiliary_labels <- make.unique(label_from_file(auxiliary, "bed"), sep = " ")
  check_unique(auxiliary_labels, "auxiliary_labels")

  check_length(coverage_labels, length(coverage), "coverage_labels", "coverage")
  check_length(coverage_colours, length(coverage), "coverage_colours", "coverage")
  check_length(coverage_groups, length(coverage), "coverage_groups", "coverage")
  if(length(coverage) > 0 && is.null(coverage_labels)) coverage_labels <- make.unique(label_from_file(coverage, "bw|bigwig"), sep = " ")
  check_unique(coverage_labels, "coverage_labels")

  missing_files <- c(xlinks, auxiliary, coverage)
  missing_files <- missing_files[!file.exists(missing_files)]
  if(length(missing_files) > 0) clipplotr_error("These input files do not exist: ", paste(missing_files, collapse = ", "))

  if(all(grepl("\\.bedgraph(\\.gz)?$", xlinks, ignore.case = TRUE))) {
    import_xlinks <- import_imaps_bedgraph
  } else if(all(grepl("\\.bed(\\.gz)?$", xlinks, ignore.case = TRUE))) {
    import_xlinks <- import.bed
  } else {
    clipplotr_error("Crosslink files need to be all in iCount bedgraph (.bedgraph) or all in bed (.bed) format.")
  }

  # Annotation
  message("INFO: Loading annotation GTF")
  annotation <- if(is.character(gtf)) load_annotation(gtf, cache_dir = cache_dir, no_ucsc = no_ucsc) else gtf

  # Crosslinks
  message("INFO: Loading crosslink track(s)")
  xlinks_grl <- lapply(xlinks, function(x) to_ucsc(import_xlinks(x), no_ucsc))

  xlink_info <- data.table(sample = labels,
                           colour = if(is.null(colours)) NA_character_ else colours,
                           group = if(is.null(groups)) NA_character_ else groups,
                           lib_size = sapply(xlinks_grl, function(x) sum(abs(x$score))),
                           size_factor = if(is.null(size_factors)) NA_real_ else size_factors)

  # Auxiliary
  auxiliary_grl <- NULL
  if(length(auxiliary) > 0) {
    message("INFO: Loading auxiliary track(s)")
    auxiliary_grl <- lapply(auxiliary, function(x) to_ucsc(import.bed(x), no_ucsc))
  }

  tracks <- list(xlinks = xlinks_grl,
                 xlink_info = xlink_info,
                 auxiliary = auxiliary_grl,
                 auxiliary_labels = auxiliary_labels,
                 coverage = coverage,
                 coverage_info = data.table(exp = coverage_labels,
                                            colour = if(is.null(coverage_colours)) NA_character_ else coverage_colours,
                                            group = if(is.null(coverage_groups)) NA_character_ else coverage_groups),
                 annotation = annotation,
                 no_ucsc = no_ucsc)
  class(tracks) <- "clipplotr_tracks"

  return(tracks)

}

#' @export
print.clipplotr_tracks <- function(x, ...) {

  cat("clipplotr tracks\n")
  cat("  crosslinks:", paste(x$xlink_info$sample, collapse = ", "), "\n")
  if(!is.null(x$auxiliary)) cat("  auxiliary: ", paste(x$auxiliary_labels, collapse = ", "), "\n")
  if(length(x$coverage) > 0) cat("  coverage:  ", paste(x$coverage_info$exp, collapse = ", "), "\n")
  invisible(x)

}

import_imaps_bedgraph <- function(bedgraph_file) {

  bg <- import.bedGraph(bedgraph_file)

  # Assign strand based on score
  strand(bg)[bg$score < 0] <- "-"
  strand(bg)[bg$score > 0] <- "+"

  # Convert scores to positives now that strands assigned
  bg$score <- abs(bg$score)

  return(bg)

}

# Read a bigWig for the region only, matching the chromosome naming style of the file
import_coverage_region <- function(file, region) {

  bw_chr <- match_seqname(as.character(seqnames(region)), seqlevels(BigWigFile(file)))

  if(is.na(bw_chr)) {

    message("WARNING: ", seqnames(region), " not found in ", file, "; coverage will be 0")
    gr <- GRanges()

  } else {

    gr <- import.bw(file, selection = BigWigSelection(GRanges(bw_chr, ranges(region))))
    seqlevels(gr) <- bw_chr
    seqlevels(gr) <- as.character(seqnames(region))

  }

  return(per_nucleotide(gr, region))

}
