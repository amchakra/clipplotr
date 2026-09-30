#' Load a GTF annotation, using a cache for later runs
#'
#' The first time a GTF is used, a TxDb and a table of genes are created and
#' saved to `cache_dir`. Later runs with the same GTF load these directly, which
#' is much faster than reading the GTF. The cache is keyed on the GTF file name,
#' size and modification time, so a changed GTF is rebuilt.
#'
#' @param gtf Path to a GTF file (can be gzip compressed).
#' @param cache_dir Directory for the cached annotation. Defaults to the user
#'   cache directory (`tools::R_user_dir("clipplotr", "cache")`).
#' @param no_ucsc If `FALSE` (default), chromosome names are converted to UCSC
#'   style (e.g. `chr1`). Set to `TRUE` for non-standard sequences.
#' @return A list with elements `txdb` (a TxDb) and `genes` (a GRanges of genes
#'   with `gene_id` and `gene_name`).
#' @export
load_annotation <- function(gtf, cache_dir = NULL, no_ucsc = FALSE) {

  if(!file.exists(gtf)) clipplotr_error("GTF file '", gtf, "' does not exist")

  if(is.null(cache_dir)) cache_dir <- tools::R_user_dir("clipplotr", which = "cache")
  info <- file.info(gtf)
  key <- sprintf("%s_%.0f_%.0f", sub("\\.(gtf|gff|gff2)(\\.gz)?$", "", basename(gtf), ignore.case = TRUE), info$size, as.numeric(info$mtime))
  txdb_file <- file.path(cache_dir, paste0(key, ".sqlite"))
  genes_file <- file.path(cache_dir, paste0(key, ".genes.rds"))

  if(file.exists(txdb_file) && file.exists(genes_file)) {

    message("INFO: ...loading cached annotation database from ", cache_dir)
    txdb <- loadDb(txdb_file)
    genes <- readRDS(genes_file)

  } else {

    message("INFO: ...creating annotation database for future runs")
    gtf_gr <- import(gtf, format = "gtf")
    genes <- gtf_gr[gtf_gr$type == "gene", c("gene_id", "gene_name")]
    txdb <- suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGRanges(gtf_gr)))

    saved <- tryCatch({
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
      saveDb(txdb, txdb_file)
      saveRDS(genes, genes_file)
      TRUE
    }, error = function(e) FALSE)
    if(!saved) message("WARNING: Could not write annotation cache to ", cache_dir, "; it will be rebuilt next time")

  }

  if(!no_ucsc) {
    seqlevelsStyle(txdb) <- "UCSC"
    seqlevelsStyle(genes) <- "UCSC"
  }

  return(list(txdb = txdb, genes = genes))

}

# Split a region into arrow segments covering each of the given features (genes or transcripts)
# so that the intron line shows the direction of transcription. About 15 arrows across the region.
arrow_segments <- function(features, ids, id_column, region) {

  spacing <- max(1, round(width(region)/15, -1), round(width(region)/15))
  region_tiled <- tile(rep(region, length(features)), width = spacing)

  segments <- GRangesList(lapply(seq_along(features), function(i) {
    gr <- subsetByOverlaps(region_tiled[[i]], features[i], type = "any")
    # Need to adjust for ends
    start(gr[1]) <- start(features[i])
    end(gr[length(gr)]) <- end(features[i])
    return(gr)
  }))

  names(segments) <- ids
  segments_dt <- as.data.table(segments)
  segments_dt[, group := NULL]
  setnames(segments_dt, "group_name", id_column)
  if(as.character(strand(region)) == "-") setnames(segments_dt, c("start", "end"), c("end", "start"))

  return(segments_dt)

}

annotation_theme <- function(p) {

  p + theme_minimal_vgrid() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
    labs(x = "Coordinate",
         y = "",
         fill = "")

}

# Gene models: all exons of each gene collapsed onto one row
plot_genes <- function(annotation, region) {

  sel_genes <- subsetByOverlaps(annotation$genes, region)

  if(length(sel_genes) == 0) {

    message("INFO: ...no genes in the region")
    return(annotation_theme(ggplot()))

  }

  sel_exons <- exons(annotation$txdb, filter = list(gene_id = sel_genes$gene_id), columns = "gene_id")

  # Region (for arrows)
  genes_order_dt <- data.table(gene_id = sel_genes$gene_id, group = seq_along(sel_genes))
  segments_dt <- arrow_segments(sel_genes, sel_genes$gene_id, "gene_id", region)
  segments_dt <- merge(segments_dt, genes_order_dt, by = "gene_id")

  # Exons
  exons_dt <- as.data.table(sel_exons)
  exons_dt$gene_id <- as.character(exons_dt$gene_id) # otherwise type is AsIs
  exons_dt <- merge(exons_dt, genes_order_dt, by = "gene_id")

  p <- ggplot() +
    geom_segment(data = segments_dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
    geom_rect(data = exons_dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene_id)) +
    scale_fill_tableau()

  return(annotation_theme(p))

}

# Transcript models: one row per transcript, coloured by gene, with thinner UTRs
plot_transcripts <- function(annotation, region) {

  txdb <- annotation$txdb

  # Get transcripts that overlap region and order for plotting
  sel_tx <- transcriptsByOverlaps(txdb, region, columns = c("gene_id", "tx_name"))
  sel_tx <- sel_tx[order(width(sel_tx), decreasing = TRUE)]

  if(length(sel_tx) == 0) {

    message("INFO: ...no transcripts in the region")
    return(annotation_theme(ggplot()))

  }

  rosetta_dt <- as.data.table(mcols(annotation$genes))[, list(gene_id, gene_name)]
  setkey(rosetta_dt, gene_id)

  tx_order_dt <- data.table(transcript_id = sel_tx$tx_name,
                            gene_id = sapply(sel_tx$gene_id, "[", 1))[, group := 1:.N]
  setkey(tx_order_dt, gene_id)
  tx_order_dt <- rosetta_dt[tx_order_dt]
  tx_order_dt[, gene := paste0(gene_name, " | ", gene_id)]

  # Region (for arrows)
  segments_dt <- arrow_segments(sel_tx, sel_tx$tx_name, "transcript_id", region)
  segments_dt <- merge(segments_dt, tx_order_dt, by = "transcript_id")

  # CDS, UTRs and, for transcripts without a CDS (e.g. ncRNA), exons
  # (half-height 0.25 for exons and CDS, 0.15 for UTRs)
  feature_table <- function(grl, tx, half_height) {

    grl <- grl[names(grl) %in% tx]
    dt <- as.data.table(grl)
    if(nrow(dt) == 0) return(NULL)
    dt[, group := NULL]
    setnames(dt, "group_name", "transcript_id")
    dt <- merge(dt, tx_order_dt, by = "transcript_id")
    dt[, `:=`(ymin = group - half_height, ymax = group + half_height)]
    return(dt)

  }

  cds_tx <- cdsBy(txdb, by = "tx", use.names = TRUE)
  noncoding_tx <- sel_tx$tx_name[!sel_tx$tx_name %in% names(cds_tx)]

  features_dt <- rbindlist(list(feature_table(exonsBy(txdb, by = "tx", use.names = TRUE), noncoding_tx, 0.25),
                                feature_table(cds_tx, sel_tx$tx_name, 0.25),
                                feature_table(fiveUTRsByTranscript(txdb, use.names = TRUE), sel_tx$tx_name, 0.15),
                                feature_table(threeUTRsByTranscript(txdb, use.names = TRUE), sel_tx$tx_name, 0.15)),
                           fill = TRUE)

  p <- ggplot() +
    geom_segment(data = segments_dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
    geom_rect(data = features_dt, mapping = aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene)) +
    scale_fill_tableau()

  return(annotation_theme(p))

}
