#' Plot CLIP tracks for one region
#'
#' @param tracks Tracks from [load_tracks()].
#' @param region A region: coordinates (`chr1:100:200:+`), gene id or gene
#'   name (see [resolve_region()]), or a GRanges.
#' @inheritParams xlink_signal
#' @param annotation Gene annotation panel: `transcript`, `gene` or `none`.
#' @param transcripts Which transcripts to show with `annotation = "transcript"`:
#'   `all`, `canonical` (Ensembl_canonical tag) or `mane` (MANE_Select tag).
#'   Genes without the tag show their canonical or longest transcript.
#' @param highlight Optional region to shade in all panels, as `c(start, end)`
#'   or `"start:end"`.
#' @param flip_x Reverse the x-axis, e.g. for genes on the negative strand.
#' @param scale_y Scale the crosslink y-axis of each group independently.
#' @param tidy_y_labels Optional approximate number of crosslink y-axis labels.
#' @param ratios Optional relative panel heights, in order: crosslinks,
#'   auxiliary, coverage, annotation. By default heights follow the number of
#'   groups, files and transcripts in each panel.
#' @param title Plot title. Defaults to the gene name, or the region as given.
#' @param base_size Base font size in points.
#' @return A patchwork object, which can be printed, saved with
#'   [ggplot2::ggsave()] or modified further. Its `"height_mm"` attribute is a
#'   suggested height in mm for a 210 mm wide figure.
#' @export
plot_region <- function(tracks, region,
                        normalisation = "libsize", smoothing = "rollmean", smoothing_window = 100,
                        annotation = "transcript", transcripts = "all", highlight = NULL, flip_x = FALSE, scale_y = FALSE,
                        tidy_y_labels = NULL, ratios = NULL, title = NULL, base_size = 10) {

  if(identical(annotation, "original")) clipplotr_error("the original annotation style has been removed; please use transcript, gene or none")
  check_choice(annotation, c("transcript", "gene", "none"), "annotation")
  check_choice(transcripts, c("all", "canonical", "mane"), "transcripts")

  if(!is.null(ratios)) {
    ratios <- suppressWarnings(as.numeric(ratios))
    if(length(ratios) != 4 || anyNA(ratios) || any(ratios < 0)) {
      clipplotr_error("ratios, please provide 4 values in order: xlink track, auxiliary tracks, coverage track, annotation track (comma separated). If a track is absent put 0.")
    }
  }

  if(!is.null(highlight)) {
    if(is.character(highlight)) highlight <- strsplit(highlight, ":")[[1]]
    highlight <- suppressWarnings(as.integer(highlight))
    if(length(highlight) != 2 || anyNA(highlight)) clipplotr_error("highlight needs to be given as start:end, e.g. 35754106:35856276")
  }

  region_name <- if(is.character(region)) region else NULL
  if(is.character(region)) region <- resolve_region(region, tracks$annotation, no_ucsc = tracks$no_ucsc)
  message("INFO: Region is ", region_string(region))

  if(is.null(title)) title <- plot_title(region_name, region, tracks$annotation)
  en_dash <- intToUtf8(0x2013)
  subtitle <- paste0(format(start(region), big.mark = ","), en_dash, format(end(region), big.mark = ","),
                     " on ", seqnames(region), " (", strand(region), " strand)")

  # Put the x-axis on the region of interest (reversed if flipped). coord_cartesian is used rather
  # than scale limits so features overlapping the region edges are kept; clip = "off" lets direct
  # labels sit in the right margin
  region_xlim <- c(start(region), end(region))
  add_x_scale <- function(p, clip = "on") {

    # At most 5 breaks as genomic coordinates make long labels
    breaks <- function(limits) {
      b <- scales::breaks_extended(n = 5)(limits)
      if(length(b) > 5) b <- b[seq(1, length(b), by = 2)]
      return(b)
    }
    labels <- scales::label_comma()
    if(flip_x) {
      p + scale_x_reverse(breaks = breaks, labels = labels) + coord_cartesian(xlim = rev(region_xlim), clip = clip)
    } else {
      p + scale_x_continuous(breaks = breaks, labels = labels) + coord_cartesian(xlim = region_xlim, clip = clip)
    }

  }

  # Shade the highlighted region underneath the data in every panel
  add_highlight <- function(p) {

    if(is.null(highlight)) return(p)
    shade <- annotate("rect", xmin = max(min(highlight), start(region)), xmax = min(max(highlight), end(region)),
                      ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.15)
    p$layers <- c(shade, p$layers)
    return(p)

  }

  message("INFO: Plotting crosslink track(s)")
  p_xlinks <- plot_xlinks(tracks, region, normalisation, smoothing, smoothing_window, scale_y, tidy_y_labels, base_size)
  panels <- list(add_x_scale(add_highlight(p_xlinks), clip = "off"))
  # Panel heights in units of 40 mm (see height_mm below)
  heights <- 0.6 * attr(p_xlinks, "n_facets")

  if(!is.null(tracks$auxiliary)) {
    message("INFO: Plotting auxiliary track(s)")
    panels <- c(panels, list(add_x_scale(add_highlight(plot_auxiliary(tracks, region, base_size)))))
    heights <- c(heights, 0.075 * length(tracks$auxiliary) + 0.05)
  }

  if(length(tracks$coverage) > 0) {
    message("INFO: Plotting coverage track(s)")
    p_coverage <- plot_coverage(tracks, region, base_size)
    panels <- c(panels, list(add_x_scale(add_highlight(p_coverage), clip = "off")))
    heights <- c(heights, 0.55 * attr(p_coverage, "n_facets"))
  }

  if(annotation != "none") {
    message("INFO: Plotting annotation tracks")
    p_annot <- if(annotation == "gene") plot_genes(tracks$annotation, region, base_size) else plot_transcripts(tracks$annotation, region, transcripts, base_size)
    panels <- c(panels, list(add_x_scale(add_highlight(p_annot))))
    heights <- c(heights, 0.12 * attr(p_annot, "n_rows") + 0.15)
  }

  if(!is.null(ratios)) {
    present <- c(TRUE, !is.null(tracks$auxiliary), length(tracks$coverage) > 0, annotation != "none")
    heights <- ratios[present]
  }

  # One coordinate axis at the bottom
  last <- length(panels)
  panels[-last] <- lapply(panels[-last], hide_x_axis)
  panels[[last]] <- panels[[last]] + labs(x = paste0(seqnames(region), " position"))

  p <- Reduce(`/`, panels) +
    plot_layout(heights = heights) +
    patchwork::plot_annotation(title = title, subtitle = subtitle, theme = theme_clipplotr(base_size))

  # Suggested height: about 40 mm per unit of panel height plus title and axis
  attr(p, "height_mm") <- round(min(400, max(100, 30 + 40 * sum(heights))))

  return(p)

}

# Gene name(s) for the title when the region was given as a gene or overlaps genes
plot_title <- function(region_name, region, annotation) {

  if(!is.null(region_name) && !grepl(":", region_name)) {
    gene <- subsetByOverlaps(annotation$genes, region, type = "equal")
    if(length(gene) == 1 && !is.na(gene$gene_name)) return(gene$gene_name)
    return(region_name)
  }

  genes <- subsetByOverlaps(annotation$genes, region)$gene_name
  genes <- unique(genes[!is.na(genes)])
  if(length(genes) == 0) return(region_string(region))
  if(length(genes) > 3) genes <- c(genes[1:3], "...")
  return(paste(genes, collapse = ", "))

}

plot_xlinks <- function(tracks, region, normalisation, smoothing, smoothing_window, scale_y, tidy_y_labels, base_size) {

  xlinks_dt <- xlink_signal(tracks, region, normalisation, smoothing, smoothing_window)
  info <- tracks$xlink_info
  colours <- if(anyNA(info$colour)) track_colours(nrow(info)) else info$colour
  grouped <- "grp" %in% colnames(xlinks_dt)

  p <- ggplot(xlinks_dt) +
    geom_line(aes(x = start, y = smoothed, group = sample, colour = sample), linewidth = 0.45) +
    scale_colour_manual(values = setNames(colours, info$sample)) +
    labs(x = "", y = signal_label(normalisation)) +
    theme_clipplotr(base_size)

  # Each group is a separate track, with its name above it
  if(grouped) p <- p + facet_wrap(~ grp, ncol = 1, scales = ifelse(scale_y, "free_y", "fixed"))
  if(!is.null(tidy_y_labels)) p <- p + scale_y_continuous(n.breaks = tidy_y_labels)

  p <- direct_labels(p, xlinks_dt, "sample", if(grouped) "grp" else NULL, base_size)
  attr(p, "n_facets") <- if(grouped) length(unique(xlinks_dt$grp)) else 1.5

  return(p)

}

plot_auxiliary <- function(tracks, region, base_size) {

  labels <- tracks$auxiliary_labels
  auxiliary_grl <- lapply(tracks$auxiliary, function(x) subsetByOverlaps(x, region, ignore.strand = FALSE))

  auxiliary_dt <- rbindlist(lapply(auxiliary_grl, as.data.table), fill = TRUE)
  auxiliary_dt$exp <- rep(labels, elementNROWS(auxiliary_grl))

  if(nrow(auxiliary_dt) == 0) message("INFO: ...no auxiliary features in the region")

  # Features without a colour (e.g. from a BED6 file mixed with BED9) are drawn in dark grey
  use_itemRgb <- "itemRgb" %in% colnames(auxiliary_dt) && nrow(auxiliary_dt) > 0
  if(use_itemRgb) {
    auxiliary_dt[is.na(itemRgb), itemRgb := "grey30"]
  } else {
    auxiliary_dt[, itemRgb := rep("grey30", .N)]
  }

  # Files in the order given, top to bottom
  auxiliary_dt[, exp := factor(exp, levels = rev(labels))]

  p <- ggplot(auxiliary_dt) +
    geom_rect(aes(xmin = start, xmax = end, ymin = as.integer(exp) - 0.35, ymax = as.integer(exp) + 0.35, fill = itemRgb)) +
    scale_fill_identity() +
    scale_y_continuous(breaks = seq_along(labels), labels = rev(labels), limits = c(0.5, length(labels) + 0.5), expand = c(0, 0)) +
    labs(x = "", y = "") +
    theme_clipplotr(base_size) +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())

  return(p)

}

plot_coverage <- function(tracks, region, base_size) {

  info <- tracks$coverage_info
  coverage_grl <- lapply(tracks$coverage, import_coverage_region, region = region)

  coverage_dt <- rbindlist(lapply(coverage_grl, as.data.table))
  coverage_dt$exp <- factor(rep(info$exp, elementNROWS(coverage_grl)), levels = info$exp)
  grouped <- !anyNA(info$group)

  if(grouped) {
    coverage_groups_dt <- data.table(exp = factor(info$exp, levels = info$exp),
                                     grp = factor(info$group, levels = unique(info$group)))
    coverage_dt <- merge(coverage_dt, coverage_groups_dt, by = "exp")
  }

  colours <- if(anyNA(info$colour)) track_colours(nrow(info)) else info$colour

  p <- ggplot(coverage_dt, aes(x = start, y = score, colour = exp)) +
    geom_line(linewidth = 0.45) +
    scale_colour_manual(values = setNames(colours, info$exp)) +
    labs(x = "", y = "Coverage") +
    theme_clipplotr(base_size)

  if(grouped) p <- p + facet_wrap(~ grp, ncol = 1, scales = "free_y")

  p <- direct_labels(p, coverage_dt, "exp", if(grouped) "grp" else NULL, base_size)
  attr(p, "n_facets") <- if(grouped) length(unique(coverage_dt$grp)) else 1.5

  return(p)

}
