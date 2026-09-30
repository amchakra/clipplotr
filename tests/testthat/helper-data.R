# The CD55 test data live in test/ at the top of the repository, which is not
# part of the built package, so tests that need them are skipped under R CMD check
test_data_dir <- function() {

  dir <- normalizePath(test_path("..", "..", "test"), mustWork = FALSE)
  skip_if_not(dir.exists(dir), "CD55 test data (test/) not available")
  return(dir)

}

test_data <- function(...) file.path(test_data_dir(), ...)

xlink_files <- function() {

  test_data(c("test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz",
              "test_hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz",
              "test_U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz",
              "test_U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz"))

}

# Annotation is cached in a temporary directory shared by all tests in a run
test_annotation <- function() {

  if(is.null(.test_env$annotation)) {
    .test_env$annotation <- suppressMessages(load_annotation(test_data("CD55_gencode.v34lift37.annotation.gtf.gz"),
                                                             cache_dir = file.path(tempdir(), "clipplotr-cache")))
  }
  return(.test_env$annotation)

}

test_tracks <- function(...) {

  suppressMessages(load_tracks(xlink_files(), gtf = test_annotation(),
                               labels = c("hnRNPC 1", "hnRNPC 2", "U2AF2 1", "U2AF2 2"),
                               groups = c("hnRNPC", "hnRNPC", "U2AF2", "U2AF2"),
                               size_factors = c(4.869687, 9.488133, 1.781117, 10.135903), ...))

}

.test_env <- new.env()
