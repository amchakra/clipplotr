#' Load a GTF annotation, using a cache for later runs
#'
#' The first time a GTF is used, a TxDb and tables of genes and transcripts are
#' created and saved to `cache_dir`. Later runs with the same GTF load these
#' directly, which is much faster than reading the GTF. The cache is keyed on
#' the GTF file name, size and modification time, so a changed GTF is rebuilt.
#'
#' @param gtf Path to a GTF file (can be gzip compressed).
#' @param cache_dir Directory for the cached annotation. Defaults to the user
#'   cache directory (`tools::R_user_dir("clipplotr", "cache")`).
#' @param no_ucsc If `FALSE` (default), chromosome names are converted to UCSC
#'   style (e.g. `chr1`). Set to `TRUE` for non-standard sequences.
#' @return A list with elements `txdb` (a TxDb), `genes` (a GRanges of genes
#'   with `gene_id` and `gene_name`) and `transcripts` (a data.table with
#'   `transcript_id`, `transcript_name`, `gene_id`, `transcript_type` and
#'   logical `canonical` and `mane` columns from the GTF tags
#'   `Ensembl_canonical` and `MANE_Select`).
#' @export
load_annotation <- function(gtf, cache_dir = NULL, no_ucsc = FALSE) {

  if(!file.exists(gtf)) clipplotr_error("GTF file '", gtf, "' does not exist")

  if(is.null(cache_dir)) cache_dir <- tools::R_user_dir("clipplotr", which = "cache")
  info <- file.info(gtf)
  # v2: cache also holds the transcript table
  key <- sprintf("%s_%.0f_%.0f_v2", sub("\\.(gtf|gff|gff2)(\\.gz)?$", "", basename(gtf), ignore.case = TRUE), info$size, as.numeric(info$mtime))
  txdb_file <- file.path(cache_dir, paste0(key, ".sqlite"))
  tables_file <- file.path(cache_dir, paste0(key, ".rds"))

  if(file.exists(txdb_file) && file.exists(tables_file)) {

    message("INFO: ...loading cached annotation database from ", cache_dir)
    txdb <- loadDb(txdb_file)
    tables <- readRDS(tables_file)

  } else {

    message("INFO: ...creating annotation database for future runs")
    gtf_gr <- import(gtf, format = "gtf")
    tables <- list(genes = gtf_gr[gtf_gr$type == "gene", c("gene_id", "gene_name")],
                   transcripts = transcript_table(gtf, gtf_gr))
    txdb <- suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGRanges(gtf_gr)))

    saved <- tryCatch({
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
      saveDb(txdb, txdb_file)
      saveRDS(tables, tables_file)
      TRUE
    }, error = function(e) FALSE)
    if(!saved) message("WARNING: Could not write annotation cache to ", cache_dir, "; it will be rebuilt next time")

  }

  genes <- tables$genes
  if(!no_ucsc) {
    seqlevelsStyle(txdb) <- "UCSC"
    seqlevelsStyle(genes) <- "UCSC"
  }

  return(list(txdb = txdb, genes = genes, transcripts = tables$transcripts))

}

# Transcript names, types and canonical/MANE flags. The tags are read from the GTF
# text because rtracklayer keeps only the last of several tag attributes on a line
transcript_table <- function(gtf, gtf_gr) {

  tx <- gtf_gr[gtf_gr$type == "transcript"]
  column <- function(name, default) if(name %in% colnames(mcols(tx))) as.character(mcols(tx)[[name]]) else default

  transcripts <- data.table(transcript_id = tx$transcript_id,
                            transcript_name = column("transcript_name", tx$transcript_id),
                            gene_id = tx$gene_id,
                            transcript_type = column("transcript_type", NA_character_))
  transcripts[is.na(transcript_name), transcript_name := transcript_id]

  lines <- readLines(gtf)
  lines <- lines[grepl("^[^\t#]*\t[^\t]*\ttranscript\t", lines)]
  tags <- data.table(transcript_id = sub('.*transcript_id "([^"]+)".*', "\\1", lines),
                     canonical = grepl('tag "Ensembl_canonical"', lines, fixed = TRUE),
                     mane = grepl('tag "MANE_Select"', lines, fixed = TRUE))
  transcripts <- merge(transcripts, unique(tags, by = "transcript_id"), by = "transcript_id", all.x = TRUE, sort = FALSE)
  transcripts[is.na(canonical), canonical := FALSE]
  transcripts[is.na(mane), mane := FALSE]

  return(transcripts)

}

# Choose which transcripts to show: all, canonical (Ensembl_canonical tag, or the
# longest transcript of genes without one) or mane (MANE_Select tag, or canonical
# for genes without one)
select_transcripts <- function(tx, transcripts, choice) {

  if(choice == "all") return(tx)

  info <- transcripts[match(tx$tx_name, transcripts$transcript_id)]
  info[, width := width(tx)]
  info[, keep := FALSE]

  if(choice == "mane") {
    info[, keep := mane]
    if(any(!info[, any(keep), by = gene_id]$V1)) message("INFO: ...genes without a MANE_Select transcript show their canonical transcript instead")
  }

  no_choice <- info[, !any(keep), by = gene_id][V1 == TRUE]$gene_id
  info[gene_id %in% no_choice, keep := canonical]

  no_canonical <- info[, !any(keep), by = gene_id][V1 == TRUE]$gene_id
  if(length(no_canonical) > 0) {
    message("INFO: ...genes without an Ensembl_canonical transcript show their longest transcript instead")
    longest <- info[gene_id %in% no_canonical, .I[which.max(width)], by = gene_id]$V1
    info[longest, keep := TRUE]
  }

  return(tx[info$keep])

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

# Shared layers for gene and transcript models: thin intron lines with chevrons
# showing the strand, boxes for exons, one labelled row per gene/transcript
gene_model_plot <- function(segments_dt, features_dt, rows_dt, base_size) {

  colours <- setNames(track_colours(length(unique(rows_dt$gene_id))), unique(rows_dt$gene_id))

  p <- ggplot() +
    geom_segment(data = segments_dt, mapping = aes(x = start, xend = end, y = group, yend = group),
                 arrow = arrow(length = unit(0.05, "in"), angle = 30, type = "open"), colour = "grey55", linewidth = 0.3) +
    geom_rect(data = features_dt, mapping = aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene_id)) +
    scale_fill_manual(values = colours) +
    scale_y_reverse(breaks = rows_dt$group, labels = rows_dt$row_label, expand = expansion(add = 0.6)) +
    labs(x = "", y = "") +
    theme_clipplotr(base_size) +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = rel(0.8)))

  attr(p, "n_rows") <- nrow(rows_dt)
  return(p)

}

empty_annotation <- function(message_text, base_size) {

  message("INFO: ...", message_text)
  p <- ggplot() + labs(x = "", y = "") + theme_clipplotr(base_size)
  attr(p, "n_rows") <- 1
  return(p)

}

# Gene models: all exons of each gene collapsed onto one row
plot_genes <- function(annotation, region, base_size = 10) {

  sel_genes <- subsetByOverlaps(annotation$genes, region)
  if(length(sel_genes) == 0) return(empty_annotation("no genes in the region", base_size))

  sel_exons <- exons(annotation$txdb, filter = list(gene_id = sel_genes$gene_id), columns = "gene_id")

  rows_dt <- data.table(gene_id = sel_genes$gene_id, group = seq_along(sel_genes),
                        row_label = ifelse(is.na(sel_genes$gene_name), sel_genes$gene_id, sel_genes$gene_name))

  segments_dt <- merge(arrow_segments(sel_genes, sel_genes$gene_id, "gene_id", region), rows_dt, by = "gene_id")

  exons_dt <- as.data.table(sel_exons)
  exons_dt$gene_id <- as.character(exons_dt$gene_id) # otherwise type is AsIs
  exons_dt <- merge(exons_dt, rows_dt, by = "gene_id")
  exons_dt[, `:=`(ymin = group - 0.25, ymax = group + 0.25)]

  return(gene_model_plot(segments_dt, exons_dt, rows_dt, base_size))

}

# Transcript models: one row per transcript, coloured by gene, with thinner UTRs
plot_transcripts <- function(annotation, region, transcripts = "all", base_size = 10) {

  txdb <- annotation$txdb

  # Get transcripts that overlap region, choose which to show and order for plotting
  sel_tx <- transcriptsByOverlaps(txdb, region, columns = c("gene_id", "tx_name"))
  if(length(sel_tx) == 0) return(empty_annotation("no transcripts in the region", base_size))
  sel_tx <- select_transcripts(sel_tx, annotation$transcripts, transcripts)

  # Rows grouped by gene (genes in order along the region), longest transcript first within each gene
  gene_ids <- sapply(sel_tx$gene_id, "[", 1)
  gene_order <- tapply(start(sel_tx), gene_ids, min)
  sel_tx <- sel_tx[order(gene_order[gene_ids], gene_ids, -width(sel_tx))]
  gene_ids <- sapply(sel_tx$gene_id, "[", 1)

  # Label each row with its transcript name, prefixed by the gene name unless the transcript
  # name already includes it (as GENCODE names do, e.g. CD55-201)
  tx_names <- annotation$transcripts$transcript_name[match(sel_tx$tx_name, annotation$transcripts$transcript_id)]
  tx_names <- ifelse(is.na(tx_names), sel_tx$tx_name, tx_names)
  gene_names <- annotation$genes$gene_name[match(gene_ids, annotation$genes$gene_id)]
  gene_names <- ifelse(is.na(gene_names), gene_ids, gene_names)
  row_labels <- ifelse(startsWith(tx_names, gene_names), tx_names, paste0(gene_names, ": ", tx_names))

  rows_dt <- data.table(transcript_id = sel_tx$tx_name,
                        gene_id = gene_ids,
                        group = seq_along(sel_tx),
                        row_label = row_labels)

  segments_dt <- merge(arrow_segments(sel_tx, sel_tx$tx_name, "transcript_id", region), rows_dt, by = "transcript_id")

  # CDS, UTRs and, for transcripts without a CDS (e.g. ncRNA), exons
  # (half-height 0.25 for exons and CDS, 0.15 for UTRs)
  feature_table <- function(grl, tx, half_height) {

    grl <- grl[names(grl) %in% tx]
    dt <- as.data.table(grl)
    if(nrow(dt) == 0) return(NULL)
    dt[, group := NULL]
    setnames(dt, "group_name", "transcript_id")
    dt <- merge(dt, rows_dt, by = "transcript_id")
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

  return(gene_model_plot(segments_dt, features_dt, rows_dt, base_size))

}
