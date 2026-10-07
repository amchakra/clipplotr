test_that("tracks are loaded with labels, groups and library sizes", {

  tracks <- test_tracks()

  expect_s3_class(tracks, "clipplotr_tracks")
  expect_equal(tracks$xlink_info$sample, c("hnRNPC 1", "hnRNPC 2", "U2AF2 1", "U2AF2 2"))
  expect_true(all(tracks$xlink_info$lib_size > 0))
  expect_equal(tracks$xlink_info$lib_size[1], sum(tracks$xlinks[[1]]$score))

})

test_that("default labels are unique file names", {

  tracks <- suppressMessages(load_tracks(xlink_files(), gtf = test_annotation()))
  expect_equal(tracks$xlink_info$sample[1], "test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events")
  expect_false(anyDuplicated(tracks$xlink_info$sample) > 0)

})

test_that("bad track options are caught before reading files", {

  annotation <- test_annotation()
  expect_error(load_tracks(xlink_files(), gtf = annotation, labels = c("a", "b")), "labels has 2 entries", class = "clipplotr_error")
  expect_error(load_tracks(xlink_files(), gtf = annotation, labels = c("a", "a", "b", "c")), "must be unique", class = "clipplotr_error")
  expect_error(load_tracks("missing.bedgraph", gtf = annotation), "do not exist", class = "clipplotr_error")
  expect_error(load_tracks(xlink_files(), gtf = annotation, size_factors = c(1, 2, 3, -1)), "positive numbers", class = "clipplotr_error")

})

test_that("coverage is read for the region only, across chromosome naming styles", {

  region <- GenomicRanges::GRanges("chr1:207513700-207513800:+")
  coverage <- import_coverage_region(test_data("test_ERR127302_plus.bigwig"), region)

  expect_equal(length(coverage), 101)
  expect_true(any(coverage$score > 0))

  # Region on a chromosome not in the bigWig gives zeros with a warning message
  expect_message(missing <- import_coverage_region(test_data("test_ERR127302_plus.bigwig"), GenomicRanges::GRanges("chr2:1-10:+")), "not found")
  expect_true(all(missing$score == 0))

})
