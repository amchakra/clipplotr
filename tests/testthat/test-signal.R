region <- GenomicRanges::GRanges("chr1:207513000-207515000:+")

test_that("libsize normalisation gives counts per million per sample", {

  tracks <- test_tracks()
  signal <- xlink_signal(tracks, region, normalisation = "libsize", smoothing = "none")

  expect_equal(levels(signal$sample), tracks$xlink_info$sample)
  expect_equal(nrow(signal), 4 * GenomicRanges::width(region))
  for(i in seq_len(nrow(tracks$xlink_info))) {
    s <- signal[sample == tracks$xlink_info$sample[i]]
    expect_equal(s$smoothed, s$score * 1e6 / tracks$xlink_info$lib_size[i])
  }

})

test_that("custom normalisation divides by size factors", {

  tracks <- test_tracks()
  signal <- xlink_signal(tracks, region, normalisation = "custom", smoothing = "none")
  s <- signal[sample == "hnRNPC 2"]
  expect_equal(s$smoothed, s$score / 9.488133)

})

test_that("rolling mean is centred and computed within each sample", {

  tracks <- test_tracks()
  raw <- xlink_signal(tracks, region, normalisation = "none", smoothing = "none")
  smoothed <- xlink_signal(tracks, region, normalisation = "none", smoothing = "rollmean", smoothing_window = 5)

  for(s in tracks$xlink_info$sample) {
    x <- raw[sample == s]$smoothed
    expected <- c(0, 0, sapply(3:(length(x) - 2), function(i) mean(x[(i - 2):(i + 2)])), 0, 0)
    expect_equal(smoothed[sample == s]$smoothed, expected)
  }

})

test_that("maxpeak scales each group to a maximum of 1, and leaves empty tracks at 0", {

  tracks <- test_tracks()
  signal <- xlink_signal(tracks, region, normalisation = "maxpeak", smoothing = "rollmean", smoothing_window = 50)
  expect_equal(signal[, max(smoothed), by = grp]$V1, c(1, 1))

  # Opposite strand has no crosslinks in the test data
  empty <- xlink_signal(tracks, GenomicRanges::GRanges("chr1:207513000-207515000:-"), normalisation = "maxpeak")
  expect_true(all(empty$smoothed == 0))

})

test_that("invalid options give clear errors", {

  tracks <- test_tracks()
  expect_error(xlink_signal(tracks, region, smoothing = "gaussian"), "gaussian smoothing has been removed", class = "clipplotr_error")
  expect_error(xlink_signal(tracks, region, normalisation = "cpm"), "normalisation needs to be one of", class = "clipplotr_error")
  expect_error(xlink_signal(tracks, GenomicRanges::GRanges("chr1:207513700-207513750:+"), smoothing_window = 100), "larger than the region", class = "clipplotr_error")

  no_size_factors <- suppressMessages(load_tracks(xlink_files(), gtf = test_annotation()))
  expect_error(xlink_signal(no_size_factors, region, normalisation = "custom"), "size factors need to be provided", class = "clipplotr_error")

})
