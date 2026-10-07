# Okabe-Ito colour-blind-safe palette, ordered so that the first colours are the
# most distinct as lines on a white background (yellow last)
okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")

# Default colours for n tracks: Okabe-Ito, or evenly spaced hues if there are more than 8
track_colours <- function(n) {

  if(n <= length(okabe_ito)) return(okabe_ito[seq_len(n)])
  return(grDevices::hcl.colors(n, palette = "Dark 3"))

}

#' The clipplotr ggplot2 theme
#'
#' A minimal theme with light horizontal grid lines, group names as headings
#' above each facet and no legend (tracks are labelled directly).
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme.
#' @export
theme_clipplotr <- function(base_size = 10) {

  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
          axis.title.y = element_text(size = rel(0.9)),
          strip.text = element_text(face = "bold", hjust = 0, size = rel(0.9)),
          strip.clip = "off",
          plot.title = element_text(face = "bold", size = rel(1.2)),
          plot.subtitle = element_text(colour = "grey40", size = rel(0.9)),
          plot.title.position = "plot",
          legend.position = "none")

}

# Remove the x-axis from all but the bottom panel so the coordinate is shown once
hide_x_axis <- function(p) {

  p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

}

# Track names written in the right margin of each facet, in the track colours,
# instead of a legend. label_column holds the names. The panel needs clip = "off"
# in its coordinates (see add_x_scale) for the labels to show.
direct_labels <- function(p, data, label_column, facet_column = NULL, base_size = 10) {

  columns <- c(label_column, facet_column)
  labels <- unique(data[, columns, with = FALSE])
  if(is.null(facet_column)) {
    labels[, position := seq_len(.N)]
  } else {
    labels[, position := seq_len(.N), by = c(facet_column)]
  }
  n_max <- max(labels$position)
  step <- min(0.25, 0.9 / n_max)
  labels[, y := 0.95 - step * (position - 1)]
  labels[, label := as.character(get(label_column))]

  # Room in the right margin for the longest label (about 0.55 em per character)
  margin_pt <- max(nchar(labels$label)) * base_size * 0.55 + 8

  p + geom_text(data = labels, mapping = aes(x = I(1.01), y = I(y), label = label, colour = .data[[label_column]]),
                hjust = 0, vjust = 1, size = base_size * 0.8 / .pt, inherit.aes = FALSE) +
    theme(plot.margin = margin(5.5, margin_pt, 5.5, 5.5))

}
