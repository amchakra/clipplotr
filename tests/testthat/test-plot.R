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

test_that("plot options are checked", {

  tracks <- test_tracks()
  expect_error(plot_region(tracks, "CD55", annotation = "original"), "original annotation style has been removed", class = "clipplotr_error")
  expect_error(plot_region(tracks, "CD55", ratios = c(1, 2)), "4 values", class = "clipplotr_error")
  expect_error(plot_region(tracks, "CD55", highlight = "100"), "start:end", class = "clipplotr_error")

})
