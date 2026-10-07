cli_options <- function() {

  list(optparse::make_option(c("-x", "--xlinks"), action = "store", type = "character", help = "Input iCLIP bedgraphs (comma separated)"),
       optparse::make_option(c("-l", "--labels"), action = "store", type = "character", help = "Unique iCLIP bedgraph labels (comma separated in same order as files)"),
       optparse::make_option(c("-c", "--colours"), action = "store", type = "character", help = "iCLIP bedgraph colours (comma separated)"),
       optparse::make_option(c("", "--groups"), action = "store", type = "character", help = "Grouping of iCLIP bedgraphs for separate plots (comma separated)"),
       optparse::make_option(c("-y", "--auxiliary"), action = "store", type = "character", help = "BED file(s) of auxiliary data (comma separated)"),
       optparse::make_option(c("", "--auxiliary_labels"), action = "store", type = "character", help = "Labels for auxiliary data (comma separated in same order as files)"),
       optparse::make_option(c("", "--coverage"), action = "store", type = "character", help = "bigwig coverage files (e.g. RNA-seq or Quantseq) - ensure same strand as region of interest (comma separated)"),
       optparse::make_option(c("", "--coverage_labels"), action = "store", type = "character", help = "Labels for coverage data (comma separated in same order as files)"),
       optparse::make_option(c("", "--coverage_colours"), action = "store", type = "character", help = "Colours for coverage data (comma separated in same order as files)"),
       optparse::make_option(c("", "--coverage_groups"), action = "store", type = "character", help = "Grouping of coverage data (comma separated in same order as files)"),
       optparse::make_option(c("", "--samples"), action = "store", type = "character", help = "Tab-separated samplesheet of tracks (columns: type, file, label, colour, group, size_factor) instead of the options above"),
       optparse::make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Reference GTF file (Gencode)"),
       optparse::make_option(c("", "--cache_dir"), action = "store", type = "character", help = "Directory for the cached annotation database [default: user cache directory]"),
       optparse::make_option(c("", "--no_ucsc"), action = "store_true", type = "logical", help = "GTF has non-standard sequences", default = FALSE),
       optparse::make_option(c("-r", "--region"), action = "store", type = "character", help = "Region(s) of interest as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b (comma separated for several)"),
       optparse::make_option(c("", "--regions"), action = "store", type = "character", help = "File of regions of interest, one per line"),
       optparse::make_option(c("", "--highlight"), action = "store", type = "character", help = "Region to highlight as 35754106:35856276"),
       optparse::make_option(c("-n", "--normalisation"), action = "store", type = "character", help = "Normalisation options: none, maxpeak, libsize, custom, libsize_maxpeak, custom_maxpeak [default %default]", default = "libsize"),
       optparse::make_option(c("", "--size_factors"), action = "store", type = "character", help = "Size factors for custom normalisation (comma separated)"),
       optparse::make_option(c("", "--scale_y"), action = "store_true", type = "logical", help = "Scale CLIP y-axis for each group independently", default = FALSE),
       optparse::make_option(c("", "--tidy_y_labels"), action = "store", type = "integer", help = "Keep this many y-axis labels (may be slightly more or fewer depending on data)"),
       optparse::make_option(c("", "--flip_x"), action = "store_true", type = "logical", help = "Flip the x-axis for plotting e.g. genes on the negative strand", default = FALSE),
       optparse::make_option(c("-s", "--smoothing"), action = "store", type = "character", help = "Smoothing options: none, rollmean [default %default]", default = "rollmean"),
       optparse::make_option(c("-w", "--smoothing_window"), action = "store", type = "integer", help = "Smoothing window [default %default]", default = 100),
       optparse::make_option(c("-a", "--annotation"), action = "store", type = "character", help = "Annotation options: gene, transcript, none [default %default]", default = "transcript"),
       optparse::make_option(c("", "--transcripts"), action = "store", type = "character", help = "Transcripts to show with transcript annotation: all, canonical (Ensembl_canonical tag) or mane (MANE_Select tag) [default %default]", default = "all"),
       optparse::make_option(c("", "--size_x"), action = "store", type = "integer", help = "Plot size in mm (x) [default: %default]", default = 210),
       optparse::make_option(c("", "--size_y"), action = "store", type = "integer", help = "Plot size in mm (y) [default: scaled to the number of tracks]"),
       optparse::make_option(c("", "--ratios"), action = "store", type = "character", help = "Specify plot ratios in order: xlink track, auxiliary tracks, coverage track, annotation track (comma separated). Put 0 if any of these track types are absent. [default: 2 for xlinks, 0.25 for 1 auxiliary track 0.5 for >1, 2 for coverage, 3 for annotation]"),
       optparse::make_option(c("-o", "--output"), action = "store", type = "character", help = "Output plot filename. For several regions include {region} in the name (one file per region) or use a .pdf (one page per region)"),
       optparse::make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE))

}

split_arg <- function(x) {

  if(is.null(x)) return(NULL)
  trimws(strsplit(x, ",")[[1]])

}

# Regions from --region (comma separated) and/or --regions (one per line, # for comments)
cli_regions <- function(opt) {

  regions <- split_arg(opt$region)
  if(!is.null(opt$regions)) {
    if(!file.exists(opt$regions)) clipplotr_error("regions file '", opt$regions, "' does not exist")
    lines <- trimws(readLines(opt$regions, warn = FALSE))
    regions <- c(regions, lines[lines != "" & !startsWith(lines, "#")])
  }
  if(length(regions) == 0) clipplotr_error("region of interest needed to be supplied as e.g. chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b")

  return(regions)

}

# File name for one region when the output contains {region}
region_output <- function(output, region) {

  gsub("{region}", gsub("[^A-Za-z0-9._+-]", "_", region), output, fixed = TRUE)

}

#' Run clipplotr from the command line
#'
#' Parses command-line arguments, loads the tracks once and writes a plot for
#' each region. Errors are printed as `ERROR: ...` and exit with status 1.
#' Run `./clipplotr --help` for all options.
#'
#' @param args Command-line arguments.
#' @return Invisibly, the output file(s) written.
#' @export
run_cli <- function(args = commandArgs(trailingOnly = TRUE)) {

  tryCatch(cli_main(args), clipplotr_error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1)
  })

}

cli_main <- function(args) {

  opt <- optparse::parse_args(optparse::OptionParser(option_list = cli_options()), args = args)
  if(opt$verbose) print(opt)

  # Checks for minimum parameters, before anything slow is loaded
  if(is.null(opt$xlinks) && is.null(opt$samples)) clipplotr_error("No iCLIP bedgraphs supplied")
  if(is.null(opt$gtf)) clipplotr_error("No Reference GTF supplied")
  if(is.null(opt$output)) clipplotr_error("No output defined")
  regions <- cli_regions(opt)

  if(length(regions) > 1 && !grepl("{region}", opt$output, fixed = TRUE) && !grepl("\\.pdf$", opt$output, ignore.case = TRUE)) {
    clipplotr_error("with several regions, --output needs to contain {region} (one file per region) or be a .pdf (one page per region)")
  }

  check_choice(opt$normalisation, normalisations, "normalisation")
  if(identical(opt$smoothing, "gaussian")) clipplotr_error("gaussian smoothing has been removed; please use rollmean or none")
  check_choice(opt$smoothing, smoothings, "smoothing")
  if(identical(opt$annotation, "original")) clipplotr_error("the original annotation style has been removed; please use transcript, gene or none")
  check_choice(opt$annotation, c("transcript", "gene", "none"), "annotation")
  check_choice(opt$transcripts, c("all", "canonical", "mane"), "transcripts")

  if(!is.null(opt$samples)) {

    track_options <- c("xlinks", "labels", "colours", "groups", "size_factors", "auxiliary", "auxiliary_labels",
                       "coverage", "coverage_labels", "coverage_colours", "coverage_groups")
    given <- track_options[track_options %in% names(opt)]
    if(length(given) > 0) clipplotr_error("--samples cannot be combined with --", paste(given, collapse = ", --"))
    track_args <- read_samplesheet(opt$samples)

  } else {

    track_args <- list(xlinks = split_arg(opt$xlinks),
                       labels = split_arg(opt$labels),
                       colours = split_arg(opt$colours),
                       groups = split_arg(opt$groups),
                       size_factors = split_arg(opt$size_factors),
                       auxiliary = split_arg(opt$auxiliary),
                       auxiliary_labels = split_arg(opt$auxiliary_labels),
                       coverage = split_arg(opt$coverage),
                       coverage_labels = split_arg(opt$coverage_labels),
                       coverage_colours = split_arg(opt$coverage_colours),
                       coverage_groups = split_arg(opt$coverage_groups))

  }

  if(opt$normalisation %in% c("custom", "custom_maxpeak") && is.null(track_args$size_factors)) {
    clipplotr_error("for custom or custom_maxpeak normalisation size factors need to be provided")
  }

  tracks <- do.call(load_tracks, c(track_args, list(gtf = opt$gtf, cache_dir = opt$cache_dir, no_ucsc = opt$no_ucsc)))

  plot_one <- function(region) {

    plot_region(tracks, region,
                normalisation = opt$normalisation,
                smoothing = opt$smoothing,
                smoothing_window = opt$smoothing_window,
                annotation = opt$annotation,
                transcripts = opt$transcripts,
                highlight = opt$highlight,
                flip_x = opt$flip_x,
                scale_y = opt$scale_y,
                tidy_y_labels = opt$tidy_y_labels,
                ratios = split_arg(opt$ratios))

  }

  # Height from --size_y, otherwise as suggested by plot_region for the tracks shown
  height <- function(p) if(is.null(opt$size_y)) attr(p, "height_mm") else opt$size_y

  if(grepl("{region}", opt$output, fixed = TRUE)) {

    outputs <- vapply(regions, region_output, character(1), output = opt$output, USE.NAMES = FALSE)
    if(anyDuplicated(outputs)) clipplotr_error("several regions give the same output file name: ", paste(unique(outputs[duplicated(outputs)]), collapse = ", "))

    for(i in seq_along(regions)) {
      if(file.exists(outputs[i])) message("WARNING: Output file '", outputs[i], "' exists and will be overwritten!")
      p <- plot_one(regions[i])
      ggsave(outputs[i], p, height = height(p), width = opt$size_x, units = "mm")
    }

  } else if(length(regions) > 1) {

    # One page per region; plots are made before opening the file so an error leaves no partial PDF
    outputs <- opt$output
    if(file.exists(outputs)) message("WARNING: Output file '", outputs, "' exists and will be overwritten!")
    plots <- lapply(regions, plot_one)
    grDevices::pdf(outputs, width = opt$size_x / 25.4, height = max(sapply(plots, height)) / 25.4)
    for(p in plots) print(p)
    grDevices::dev.off()

  } else {

    outputs <- opt$output
    if(file.exists(outputs)) message("WARNING: Output file '", outputs, "' exists and will be overwritten!")
    p <- plot_one(regions)
    ggsave(outputs, p, height = height(p), width = opt$size_x, units = "mm")

  }

  message("Completed")
  invisible(outputs)

}
