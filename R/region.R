#' Resolve a region of interest
#'
#' @param region Coordinates as `chromosome:start:end:strand` (or
#'   `chromosome:start-end:strand`), an Ensembl/GENCODE gene id (with or
#'   without version suffix), or a gene name.
#' @param annotation Annotation from [load_annotation()], used to look up gene
#'   ids and names.
#' @param no_ucsc If `FALSE` (default), chromosome names in coordinates are
#'   converted to UCSC style where possible.
#' @return A GRanges of length 1 with a strand.
#' @export
resolve_region <- function(region, annotation, no_ucsc = FALSE) {

  genes <- annotation$genes
  coords <- regmatches(region, regexec("^([^:]+):([0-9]+)[:-]([0-9]+):([+-])$", region))[[1]]

  if(length(coords) == 5) {

    if(as.numeric(coords[4]) < as.numeric(coords[3])) clipplotr_error("region '", region, "' end is before its start")
    region_gr <- GRanges(seqnames = coords[2],
                         ranges = IRanges(start = as.integer(coords[3]), end = as.integer(coords[4])),
                         strand = coords[5])

    # Non-standard contigs have no UCSC equivalent, so leave their name as given
    if(!no_ucsc) region_gr <- tryCatch(to_ucsc(region_gr, no_ucsc), error = function(e) region_gr)

  } else if(grepl(":", region)) {

    clipplotr_error("region '", region, "' should be given as chromosome:start:end:strand, e.g. chr3:35754106:35856276:+")

  } else {

    if(grepl("^ENS", region)) {

      # Match gene id with or without the version suffix
      unversioned <- function(x) sub("\\..*$", "", x)
      region_gr <- genes[genes$gene_id == region | unversioned(genes$gene_id) == unversioned(region)]

    } else {

      region_gr <- genes[genes$gene_name %in% region]

    }

    if(length(region_gr) == 0) clipplotr_error("region '", region, "' was not found in the GTF")
    if(length(region_gr) > 1) {
      clipplotr_error("There are more than one genes with the identifier '", region, "' (", paste(region_gr$gene_id, collapse = ", "), "). Consider specifying the coordinates instead.")
    }

  }

  seqlevels(region_gr) <- as.character(seqnames(region_gr)) # Cut down to one seqlevel for later comparison
  mcols(region_gr) <- NULL

  return(region_gr)

}
