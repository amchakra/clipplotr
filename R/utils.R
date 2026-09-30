# Errors raised with this class are user-facing: the command-line tool prints
# them as "ERROR: ..." and exits with status 1 rather than showing a traceback
clipplotr_error <- function(...) {

  stop(structure(class = c("clipplotr_error", "error", "condition"),
                 list(message = paste0(...), call = NULL)))

}

# Check that a per-file option has one entry per file
check_length <- function(values, n, option, files_option) {

  if(!is.null(values) && length(values) != n) {
    clipplotr_error(option, " has ", length(values), " entries but ", files_option, " has ", n, " files")
  }

}

check_unique <- function(values, option) {

  if(anyDuplicated(values)) {
    clipplotr_error(option, " must be unique, but these are duplicated: ", paste(unique(values[duplicated(values)]), collapse = ", "))
  }

}

check_choice <- function(value, choices, option) {

  if(length(value) != 1 || !value %in% choices) {
    clipplotr_error(option, " needs to be one of ", paste(choices, collapse = ", "))
  }

}

# Default track label: file name without directory or extension
label_from_file <- function(files, extensions) {

  sub(paste0("\\.(", extensions, ")(\\.gz)?$"), "", basename(files), ignore.case = TRUE)

}

# Find the name used for a chromosome in a file with a different naming style (e.g. chr1 vs 1)
match_seqname <- function(chr, available) {

  candidates <- unique(c(chr, sub("^chr", "", chr), paste0("chr", chr)))
  if(chr %in% c("chrM", "M", "MT")) candidates <- c(candidates, "chrM", "MT")
  hit <- candidates[candidates %in% available]
  if(length(hit) == 0) return(NA_character_)
  return(hit[1])

}

# Convert to UCSC chromosome names unless asked not to
to_ucsc <- function(x, no_ucsc) {

  if(!no_ucsc) seqlevelsStyle(x) <- "UCSC"
  return(x)

}

# One score per nucleotide across the region, 0 where the track has no signal
per_nucleotide <- function(gr, region) {

  positions <- unlist(tile(region, width = 1))

  ol <- findOverlaps(positions, gr)
  positions$score <- 0
  positions[queryHits(ol)]$score <- gr[subjectHits(ol)]$score

  return(positions)

}

region_string <- function(region) {

  paste0(seqnames(region), ":", start(region), "-", end(region), ":", strand(region))

}
