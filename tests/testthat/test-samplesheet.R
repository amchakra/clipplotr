write_sheet <- function(lines) {

  file <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
  writeLines(lines, file)
  return(file)

}

test_that("samplesheet rows become load_tracks arguments", {

  file <- write_sheet(c("type\tfile\tlabel\tcolour\tgroup\tsize_factor",
                        "xlink\ta.bedgraph.gz\tA 1\t#000000\tA\t1.5",
                        "xlink\t/abs/b.bedgraph.gz\tA 2\t#111111\tA\t2",
                        "auxiliary\tpeaks.bed.gz\tpeaks\t\t\t",
                        "coverage\tc.bigwig\tRNA\t\t\t"))
  args <- read_samplesheet(file)

  expect_equal(args$xlinks, c(file.path(dirname(file), "a.bedgraph.gz"), "/abs/b.bedgraph.gz"))
  expect_equal(args$labels, c("A 1", "A 2"))
  expect_equal(args$colours, c("#000000", "#111111"))
  expect_equal(args$groups, c("A", "A"))
  expect_equal(args$size_factors, c("1.5", "2"))
  expect_equal(args$auxiliary_labels, "peaks")
  expect_equal(args$coverage, file.path(dirname(file), "c.bigwig"))
  expect_null(args$coverage_colours)

})

test_that("optional columns can be left out", {

  args <- read_samplesheet(write_sheet(c("type\tfile", "xlink\ta.bedgraph", "xlink\tb.bedgraph")))
  expect_length(args$xlinks, 2)
  expect_null(args$labels)
  expect_null(args$auxiliary)

})

test_that("bad samplesheets give clear errors", {

  expect_error(read_samplesheet(write_sheet(c("file", "a.bedgraph"))), "missing column", class = "clipplotr_error")
  expect_error(read_samplesheet(write_sheet(c("type\tfile\tcolor", "xlink\ta.bedgraph\tred"))), "unknown column", class = "clipplotr_error")
  expect_error(read_samplesheet(write_sheet(c("type\tfile", "peak\ta.bed"))), "xlink, auxiliary or coverage", class = "clipplotr_error")
  expect_error(read_samplesheet(write_sheet(c("type\tfile", "coverage\ta.bw"))), "at least one xlink", class = "clipplotr_error")
  expect_error(read_samplesheet(write_sheet(c("type\tfile\tlabel", "xlink\ta.bedgraph\tA", "xlink\tb.bedgraph\t"))), "filled for all xlink rows", class = "clipplotr_error")

})
