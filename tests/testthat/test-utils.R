test_that("default labels drop directory and extension", {

  expect_equal(label_from_file(c("a/b/x.bedgraph.gz", "y.bedGraph", "z.bed"), "bedgraph|bed"), c("x", "y", "z"))
  expect_equal(label_from_file("dir/sample.bw", "bw|bigwig"), "sample")

})

test_that("chromosome names are matched across naming styles", {

  expect_equal(match_seqname("chr1", c("1", "2")), "1")
  expect_equal(match_seqname("1", c("chr1", "chr2")), "chr1")
  expect_equal(match_seqname("chrM", c("MT")), "MT")
  expect_true(is.na(match_seqname("chr1", c("chr2"))))

})

test_that("per-nucleotide scores fill gaps with 0", {

  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 5), strand = "+")
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2, 4), width = 1), strand = "+", score = c(3, 7))
  expect_equal(per_nucleotide(gr, region)$score, c(0, 3, 0, 7, 0))

})

test_that("per-file option checks give clear errors", {

  expect_error(check_length(c("a", "b"), 3, "labels", "xlinks"), "labels has 2 entries but xlinks has 3 files", class = "clipplotr_error")
  expect_error(check_unique(c("a", "a", "b"), "labels"), "duplicated: a", class = "clipplotr_error")
  expect_silent(check_length(NULL, 3, "labels", "xlinks"))

})
