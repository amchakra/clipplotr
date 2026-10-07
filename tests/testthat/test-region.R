test_that("coordinates are parsed in both formats", {

  annotation <- test_annotation()

  region <- resolve_region("chr1:207513000:207515000:+", annotation)
  expect_equal(c(start(region), end(region)), c(207513000, 207515000))
  expect_equal(as.character(strand(region)), "+")

  expect_equal(resolve_region("chr1:207513000-207515000:-", annotation), GenomicRanges::GRanges("chr1:207513000-207515000:-"))

})

test_that("genes are found by name and by id with or without version", {

  annotation <- test_annotation()
  by_name <- resolve_region("CD55", annotation)

  expect_equal(resolve_region("ENSG00000196352.16_8", annotation), by_name)
  expect_equal(resolve_region("ENSG00000196352", annotation), by_name)
  expect_equal(as.character(strand(by_name)), "+")

})

test_that("bad regions give clear errors", {

  annotation <- test_annotation()

  expect_error(resolve_region("NOTAGENE", annotation), "was not found", class = "clipplotr_error")
  expect_error(resolve_region("ENSG0000019635", annotation), "was not found", class = "clipplotr_error")
  expect_error(resolve_region("chr1:100:200", annotation), "should be given as", class = "clipplotr_error")
  expect_error(resolve_region("chr1:200:100:+", annotation), "end is before its start", class = "clipplotr_error")

})
