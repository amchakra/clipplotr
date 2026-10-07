test_that("plot_region returns one panel per track type", {

  tracks <- test_tracks(auxiliary = test_data("test_Alu_rev.bed.gz"),
                        coverage = test_data(c("test_ERR127302_plus.bigwig", "test_ERR127306_plus.bigwig")))

  p <- suppressMessages(plot_region(tracks, "CD55", annotation = "transcript", highlight = "207513650:207513800"))
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 3) # auxiliary, coverage and annotation below the crosslinks

  p <- suppressMessages(plot_region(tracks, "chr1:207513000:207515000:+", annotation = "none", flip_x = TRUE))
  expect_length(p$patches$plots, 2)

  # Rendering works for each annotation style
  for(annotation in c("transcript", "gene", "none")) {
    p <- suppressMessages(plot_region(tracks, "chr1:207513000:207515000:+", annotation = annotation, smoothing_window = 50))
    file <- withr::local_tempfile(fileext = ".png")
    expect_no_error(suppressMessages(ggplot2::ggsave(file, p, width = 200, height = 250, units = "mm")))
    expect_gt(file.size(file), 0)
  }

})

test_that("title, labels, colours and suggested height follow the tracks", {

  tracks <- test_tracks()

  p <- suppressMessages(plot_region(tracks, "CD55", annotation = "none"))
  expect_equal(p$patches$annotation$title, "CD55")
  expect_match(p$patches$annotation$subtitle, "on chr1 \\(\\+ strand\\)")

  # Coordinates overlapping a gene are titled with the gene name
  p <- suppressMessages(plot_region(tracks, "chr1:207513000:207515000:+", annotation = "none"))
  expect_equal(p$patches$annotation$title, "CD55")

  # No legend: tracks are labelled in the margin with the default Okabe-Ito colours
  xlinks_panel <- p[[1]]
  labels <- ggplot2::layer_data(xlinks_panel, length(xlinks_panel$layers))
  expect_setequal(labels$label, tracks$xlink_info$sample)
  expect_setequal(labels$colour, track_colours(4))

  # More tracks make a taller figure
  short <- attr(suppressMessages(plot_region(tracks, "CD55", annotation = "none")), "height_mm")
  tall <- attr(suppressMessages(plot_region(tracks, "CD55", annotation = "transcript")), "height_mm")
  expect_gt(tall, short)

})

test_that("default colours are colour-blind safe and extend beyond 8 tracks", {

  expect_equal(track_colours(3), c("#0072B2", "#D55E00", "#009E73"))
  expect_length(unique(track_colours(12)), 12)

})

test_that("plot options are checked", {

  tracks <- test_tracks()
  expect_error(plot_region(tracks, "CD55", annotation = "original"), "original annotation style has been removed", class = "clipplotr_error")
  expect_error(plot_region(tracks, "CD55", transcripts = "longest"), "transcripts needs to be one of", class = "clipplotr_error")
  expect_error(plot_region(tracks, "CD55", ratios = c(1, 2)), "4 values", class = "clipplotr_error")
  expect_error(plot_region(tracks, "CD55", highlight = "100"), "start:end", class = "clipplotr_error")

})
