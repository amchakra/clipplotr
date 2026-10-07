normalisations <- c("none", "maxpeak", "libsize", "custom", "custom_maxpeak", "libsize_maxpeak")
smoothings <- c("none", "rollmean")

#' Normalised and smoothed crosslink signal for a region
#'
#' @param tracks Tracks from [load_tracks()].
#' @param region A GRanges from [resolve_region()].
#' @param normalisation One of `none`, `maxpeak`, `libsize` (counts per
#'   million), `custom` (divide by size factors), `libsize_maxpeak` or
#'   `custom_maxpeak`. The `maxpeak` options scale each sample (or each group,
#'   if groups are given) to its maximum after smoothing.
#' @param smoothing `rollmean` (centred rolling mean) or `none`.
#' @param smoothing_window Rolling mean window in nucleotides.
#' @return A data.table with one row per sample and nucleotide.
#' @export
xlink_signal <- function(tracks, region, normalisation = "libsize", smoothing = "rollmean", smoothing_window = 100) {

  check_choice(normalisation, normalisations, "normalisation")
  if(identical(smoothing, "gaussian")) clipplotr_error("gaussian smoothing has been removed; please use rollmean or none")
  check_choice(smoothing, smoothings, "smoothing")
  if(!is.numeric(smoothing_window) || smoothing_window < 1) clipplotr_error("smoothing_window needs to be a positive integer")
  if(smoothing == "rollmean" && smoothing_window > width(region)) {
    clipplotr_error("smoothing_window (", smoothing_window, ") is larger than the region (", width(region), " nt)")
  }

  info <- tracks$xlink_info
  if(normalisation %in% c("custom", "custom_maxpeak")) {
    if(anyNA(info$size_factor)) clipplotr_error("for custom or custom_maxpeak normalisation size factors need to be provided")
    info$libSize <- info$size_factor
  } else {
    info$libSize <- info$lib_size
  }

  # Subset for region and add in 0 count positions
  xlinks_dt <- rbindlist(lapply(seq_along(tracks$xlinks), function(i) {

    dt <- as.data.table(per_nucleotide(tracks$xlinks[[i]], region))
    dt$sample <- info$sample[i]
    dt$libSize <- info$libSize[i]
    return(dt)

  }))
  xlinks_dt[, sample := factor(sample, levels = info$sample)]

  # Do the normalisation
  xlinks_dt[, norm := switch(normalisation,
                             "libsize" = (score * 1e6)/libSize,
                             "libsize_maxpeak" = (score * 1e6)/libSize,
                             "maxpeak" = score,
                             "none" = score,
                             "custom" = score/libSize,
                             "custom_maxpeak" = score/libSize),
            by = sample]

  # Do the smoothing
  xlinks_dt[, smoothed := switch(smoothing,
                                 "rollmean" = frollmean(norm, smoothing_window, align = "center", fill = 0),
                                 "none" = norm),
            by = sample]

  # Assign groups
  if(!anyNA(info$group)) {

    groups_dt <- data.table(sample = factor(info$sample, levels = info$sample),
                            grp = factor(info$group, levels = unique(info$group)))
    xlinks_dt <- merge(xlinks_dt, groups_dt, by = "sample")

  }

  # Now normalise to max peak after smoothing (leave tracks with no signal at 0)
  scale_to_max <- function(x) if(max(x) > 0) x/max(x) else x
  if(grepl("maxpeak", normalisation)) {

    if("grp" %in% colnames(xlinks_dt)) {
      xlinks_dt[, smoothed := scale_to_max(smoothed), by = grp]
    } else {
      xlinks_dt[, smoothed := scale_to_max(smoothed), by = sample]
    }

  }

  return(xlinks_dt[])

}

signal_label <- function(normalisation) {

  switch(normalisation,
         "libsize" = "Crosslink signal\n(counts per million)",
         "libsize_maxpeak" = "Crosslink signal\n(counts per million\nnormalised to max peak)",
         "maxpeak" = "Crosslink signal\n(normalised to max peak)",
         "none" = "Raw crosslink signal",
         "custom" = "Crosslink signal\n(custom normalisation)",
         "custom_maxpeak" = "Crosslink signal\n(custom normalisation\nadjusted to max peak)")

}
