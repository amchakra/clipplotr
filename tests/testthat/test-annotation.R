# Two genes: G1 has a MANE and a canonical-only transcript plus a longer untagged one,
# G2 has no tags at all
tagged_gtf <- function() {

  file <- withr::local_tempfile(fileext = ".gtf", .local_envir = parent.frame())
  writeLines(c(
    'chr1\tX\tgene\t100\t900\t.\t+\t.\tgene_id "G1"; gene_name "A";',
    'chr1\tX\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_name "A-201"; transcript_type "protein_coding"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select";',
    'chr1\tX\texon\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number 1;',
    'chr1\tX\ttranscript\t100\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; transcript_name "A-202"; transcript_type "protein_coding"; tag "basic";',
    'chr1\tX\texon\t100\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; exon_number 1;',
    'chr1\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T3"; transcript_name "A-203"; transcript_type "retained_intron";',
    'chr1\tX\texon\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T3"; exon_number 1;',
    'chr1\tX\tgene\t1000\t1400\t.\t+\t.\tgene_id "G2"; gene_name "B";',
    'chr1\tX\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "G2"; transcript_id "T4"; transcript_name "B-201";',
    'chr1\tX\texon\t1000\t1200\t.\t+\t.\tgene_id "G2"; transcript_id "T4"; exon_number 1;',
    'chr1\tX\ttranscript\t1000\t1400\t.\t+\t.\tgene_id "G2"; transcript_id "T5"; transcript_name "B-202";',
    'chr1\tX\texon\t1000\t1400\t.\t+\t.\tgene_id "G2"; transcript_id "T5"; exon_number 1;'
  ), file)
  return(file)

}

test_that("transcript tags are read even when a line has several", {

  annotation <- suppressMessages(load_annotation(tagged_gtf(), cache_dir = withr::local_tempdir()))
  tx <- annotation$transcripts[order(transcript_id)]

  expect_equal(tx$transcript_name, c("A-201", "A-202", "A-203", "B-201", "B-202"))
  expect_equal(tx$canonical, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(tx$mane, c(TRUE, FALSE, FALSE, FALSE, FALSE))

})

test_that("canonical and MANE selection fall back to the longest transcript", {

  annotation <- suppressMessages(load_annotation(tagged_gtf(), cache_dir = withr::local_tempdir()))
  region <- GenomicRanges::GRanges("chr1:1-2000:+")
  tx <- GenomicFeatures::transcriptsByOverlaps(annotation$txdb, region, columns = c("gene_id", "tx_name"))

  expect_setequal(select_transcripts(tx, annotation$transcripts, "all")$tx_name, paste0("T", 1:5))
  expect_message(canonical <- select_transcripts(tx, annotation$transcripts, "canonical"), "longest transcript")
  expect_setequal(canonical$tx_name, c("T1", "T5"))
  expect_setequal(suppressMessages(select_transcripts(tx, annotation$transcripts, "mane"))$tx_name, c("T1", "T5"))

})

test_that("transcript rows are labelled with transcript names", {

  annotation <- suppressMessages(load_annotation(tagged_gtf(), cache_dir = withr::local_tempdir()))
  region <- GenomicRanges::GRanges("chr1:1-2000:+")

  p <- plot_transcripts(annotation, region)
  expect_equal(attr(p, "n_rows"), 5)
  # Grouped by gene, longest first within each gene
  expect_equal(ggplot2::layer_scales(p)$y$labels, c("A-203", "A-202", "A-201", "B-202", "B-201"))

  p <- suppressMessages(plot_transcripts(annotation, region, transcripts = "canonical"))
  expect_equal(attr(p, "n_rows"), 2)

  p <- plot_genes(annotation, region)
  expect_equal(ggplot2::layer_scales(p)$y$labels, c("A", "B"))

})

test_that("transcripts without a name that includes the gene are prefixed with the gene name", {

  file <- withr::local_tempfile(fileext = ".gtf")
  writeLines(c(
    'chr1\tX\tgene\t100\t900\t.\t+\t.\tgene_id "G1"; gene_name "A";',
    'chr1\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
    'chr1\tX\texon\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
    'chr1\tX\tgene\t50\t400\t.\t+\t.\tgene_id "G2"; gene_name "B";',
    'chr1\tX\ttranscript\t50\t400\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; transcript_name "isoform1";',
    'chr1\tX\texon\t50\t400\t.\t+\t.\tgene_id "G2"; transcript_id "T2";'
  ), file)
  annotation <- suppressMessages(load_annotation(file, cache_dir = withr::local_tempdir()))

  # B starts first along the region, so its rows come first
  p <- plot_transcripts(annotation, GenomicRanges::GRanges("chr1:1-1000:+"))
  expect_equal(ggplot2::layer_scales(p)$y$labels, c("B: isoform1", "A: T1"))

})
