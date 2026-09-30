#' Read a samplesheet of tracks
#'
#' A tab-separated file with one row per track, as an alternative to passing
#' comma-separated lists of files, labels, colours and groups. Columns:
#'
#' * `type` (required): `xlink`, `auxiliary` or `coverage`
#' * `file` (required): path to the file. Relative paths are relative to the
#'   samplesheet's directory.
#' * `label`, `colour`, `group`, `size_factor` (optional): as for the
#'   corresponding [load_tracks()] arguments. `colour` and `group` apply to
#'   `xlink` and `coverage` rows, `size_factor` to `xlink` rows. Leave a cell
#'   empty to use the default.
#'
#' @param file Path to the samplesheet.
#' @return A named list of arguments for [load_tracks()] (without `gtf`).
#' @export
read_samplesheet <- function(file) {

  if(!file.exists(file)) clipplotr_error("samplesheet '", file, "' does not exist")

  sheet <- fread(file, sep = "\t", colClasses = "character", na.strings = "", fill = TRUE)
  setnames(sheet, tolower(trimws(colnames(sheet))))

  missing_columns <- setdiff(c("type", "file"), colnames(sheet))
  if(length(missing_columns) > 0) clipplotr_error("samplesheet is missing column(s): ", paste(missing_columns, collapse = ", "))

  unknown_columns <- setdiff(colnames(sheet), c("type", "file", "label", "colour", "group", "size_factor"))
  if(length(unknown_columns) > 0) clipplotr_error("samplesheet has unknown column(s): ", paste(unknown_columns, collapse = ", "))

  sheet[, type := tolower(trimws(type))]
  unknown_types <- setdiff(sheet$type, c("xlink", "auxiliary", "coverage"))
  if(length(unknown_types) > 0) clipplotr_error("samplesheet type needs to be xlink, auxiliary or coverage, not: ", paste(unknown_types, collapse = ", "))
  if(anyNA(sheet$file)) clipplotr_error("samplesheet has rows without a file")

  # Relative paths are relative to the samplesheet
  relative <- !grepl("^(/|~|[A-Za-z]:)", sheet$file)
  sheet$file[relative] <- file.path(dirname(file), sheet$file[relative])

  # A column for rows of one type: NULL if absent or entirely empty, otherwise every row must be filled
  column <- function(rows, name) {

    if(!name %in% colnames(sheet)) return(NULL)
    values <- sheet[[name]][rows]
    if(all(is.na(values))) return(NULL)
    if(anyNA(values)) clipplotr_error("samplesheet column '", name, "' needs to be filled for all ", unique(sheet$type[rows]), " rows or none")
    return(values)

  }

  xlink <- sheet$type == "xlink"
  auxiliary <- sheet$type == "auxiliary"
  coverage <- sheet$type == "coverage"
  if(!any(xlink)) clipplotr_error("samplesheet needs at least one xlink row")

  args <- list(xlinks = sheet$file[xlink],
               labels = column(xlink, "label"),
               colours = column(xlink, "colour"),
               groups = column(xlink, "group"),
               size_factors = column(xlink, "size_factor"),
               auxiliary = if(any(auxiliary)) sheet$file[auxiliary] else NULL,
               auxiliary_labels = column(auxiliary, "label"),
               coverage = if(any(coverage)) sheet$file[coverage] else NULL,
               coverage_labels = column(coverage, "label"),
               coverage_colours = column(coverage, "colour"),
               coverage_groups = column(coverage, "group"))

  return(args)

}
