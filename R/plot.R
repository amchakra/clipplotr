#' Plot CLIP tracks for one region
#'
#' @param tracks Tracks from [load_tracks()].
#' @param region A region: coordinates (`chr1:100:200:+`), gene id or gene
#'   name (see [resolve_region()]), or a GRanges.
#' @inheritParams xlink_signal
#' @param annotation Gene annotation panel: `transcript`, `gene` or `none`.
#' @param highlight Optional region to shade, as `c(start, end)` or `"start:end"`.
#' @param flip_x Reverse the x-axis, e.g. for genes on the negative strand.
#' @param scale_y Scale the crosslink y-axis of each group independently.
#' @param tidy_y_labels Optional approximate number of crosslink y-axis labels.
#' @param ratios Optional relative panel heights, in order: crosslinks,
#'   auxiliary, coverage, annotation. Defaults to 2, 0.25 (0.5 for more than one
#'   auxiliary file), 2 and 3.
#' @param title Plot title. Defaults to the region as given.
#' @return A patchwork object, which can be printed, saved with
#'   [ggplot2::ggsave()] or modified further.
#' @export
plot_region <- function(tracks, region,
                        normalisation = "libsize", smoothing = "rollmean", smoothing_window = 100,
                        annotation = "transcript", highlight = NULL, flip_x = FALSE, scale_y = FALSE,
                        tidy_y_labels = NULL, ratios = NULL, title = NULL) {

  if(identical(annotation, "original")) clipplotr_error("the original annotation style has been removed; please use transcript, gene or none")
  check_choice(annotation, c("transcript", "gene", "none"), "annotation")

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

  if(is.null(title)) title <- if(is.character(region)) region else region_string(region)
  if(is.character(region)) region <- resolve_region(region, tracks$annotation, no_ucsc = tracks$no_ucsc)
  message("INFO: Region is ", region_string(region))

  # Put the x-axis on the region of interest (reversed if flipped). coord_cartesian is used
  # rather than scale limits so that features overlapping the region edges are kept
  add_x_scale <- function(p) {

    region_xlim <- c(start(region), end(region))
    if(flip_x) {
      p + scale_x_reverse() + coord_cartesian(xlim = rev(region_xlim))
    } else {
      p + scale_x_continuous() + coord_cartesian(xlim = region_xlim)
    }

  }

  message("INFO: Plotting crosslink track(s)")
  panels <- list(add_x_scale(plot_xlinks(tracks, region, normalisation, smoothing, smoothing_window, highlight, scale_y, tidy_y_labels, title)))
  heights <- if(is.null(ratios)) 2 else ratios[1]

  if(!is.null(tracks$auxiliary)) {
    message("INFO: Plotting auxiliary track(s)")
    panels <- c(panels, list(add_x_scale(plot_auxiliary(tracks, region))))
    heights <- c(heights, if(is.null(ratios)) ifelse(length(tracks$auxiliary) > 1, 0.5, 0.25) else ratios[2])
  }

  if(length(tracks$coverage) > 0) {
    message("INFO: Plotting coverage track(s)")
    panels <- c(panels, list(add_x_scale(plot_coverage(tracks, region))))
    heights <- c(heights, if(is.null(ratios)) 2 else ratios[3])
  }

  if(annotation != "none") {
    message("INFO: Plotting annotation tracks")
    p_annot <- if(annotation == "gene") plot_genes(tracks$annotation, region) else plot_transcripts(tracks$annotation, region)
    panels <- c(panels, list(add_x_scale(p_annot)))
    heights <- c(heights, if(is.null(ratios)) 3 else ratios[4])
  }

  p <- Reduce(`/`, panels) + plot_layout(heights = heights)

  return(p)

}

plot_xlinks <- function(tracks, region, normalisation, smoothing, smoothing_window, highlight, scale_y, tidy_y_labels, title) {

  xlinks_dt <- xlink_signal(tracks, region, normalisation, smoothing, smoothing_window)
  info <- tracks$xlink_info

  p <- ggplot(xlinks_dt) +
    geom_line(aes(x = start, y = smoothed, group = sample, color = sample)) +
    labs(title = title,
         x = "",
         y = signal_label(normalisation),
         colour = "") +
    theme_minimal_grid() + theme(legend.position = "top")

  if(anyNA(info$colour)) {
    p <- p + scale_colour_tableau(palette = "Tableau 10")
  } else {
    p <- p + scale_colour_manual(values = setNames(info$colour, info$sample))
  }

  # Facet if groups
  if("grp" %in% colnames(xlinks_dt)) {
    p <- p +
      facet_grid(grp ~ ., scales = ifelse(scale_y, "free_y", "fixed"), labeller = label_wrap_gen(10)) +
      theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0))
  }

  # Add highlight
  if(!is.null(highlight)) {
    highlight_dt <- data.table(x1 = highlight[1],
                               x2 = highlight[2],
                               y1 = 0,
                               y2 = Inf) # Assumes not facetted with free scales
    p <- p + geom_rect(data = highlight_dt, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = "grey50", alpha = 0.2)
  }

  if(!is.null(tidy_y_labels)) p <- p + scale_y_continuous(n.breaks = tidy_y_labels)

  return(p)

}

plot_auxiliary <- function(tracks, region) {

  labels <- tracks$auxiliary_labels
  auxiliary_grl <- lapply(tracks$auxiliary, function(x) subsetByOverlaps(x, region, ignore.strand = FALSE))

  auxiliary_dt <- rbindlist(lapply(auxiliary_grl, as.data.table), fill = TRUE)
  auxiliary_dt$exp <- rep(labels, elementNROWS(auxiliary_grl))
  auxiliary_dt[, centre := start + width/2]

  if(nrow(auxiliary_dt) == 0) message("INFO: ...no auxiliary features in the region")

  # Features without a colour (e.g. from a BED6 file mixed with BED9) are drawn in the default dark grey
  use_itemRgb <- "itemRgb" %in% colnames(auxiliary_dt) && nrow(auxiliary_dt) > 0
  if(use_itemRgb) auxiliary_dt[is.na(itemRgb), itemRgb := "grey20"]

  p <- ggplot(auxiliary_dt, aes(x = centre, width = width, y = exp)) +
    scale_y_discrete(breaks = labels,
                     limits = rev(sort(labels))) +
    labs(y = "",
         x = "") +
    theme_minimal_grid()

  if(use_itemRgb) {

    auxiliary_cols <- unique(auxiliary_dt$itemRgb)
    names(auxiliary_cols) <- auxiliary_cols

    p <- p +
      geom_tile(aes(fill = itemRgb)) +
      scale_fill_manual(values = auxiliary_cols) +
      theme(legend.position = "none")

  } else {

    p <- p + geom_tile()

  }

  return(p)

}

plot_coverage <- function(tracks, region) {

  info <- tracks$coverage_info
  coverage_grl <- lapply(tracks$coverage, import_coverage_region, region = region)

  coverage_dt <- rbindlist(lapply(coverage_grl, as.data.table))
  coverage_dt$exp <- factor(rep(info$exp, elementNROWS(coverage_grl)), levels = info$exp)

  if(!anyNA(info$group)) {
    coverage_groups_dt <- data.table(exp = factor(info$exp, levels = info$exp),
                                     grp = factor(info$group, levels = unique(info$group)))
    coverage_dt <- merge(coverage_dt, coverage_groups_dt, by = "exp")
  }

  p <- ggplot(coverage_dt, aes(x = start, y = score, col = exp)) +
    geom_line() +
    labs(x = "",
         y = "Coverage",
         colour = "") +
    theme_minimal_grid() + theme(legend.position = "bottom") + theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0))

  if(!anyNA(info$colour)) {
    p <- p + scale_colour_manual(values = setNames(info$colour, info$exp))
  } else {
    p <- p + scale_colour_tableau(palette = "Superfishel Stone")
  }

  if("grp" %in% colnames(coverage_dt)) {
    p <- p + facet_grid(grp ~ ., scales = "free_y", labeller = label_wrap_gen(10)) + theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0))
  }

  return(p)

}
