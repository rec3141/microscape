#' @title Load Sample Metadata
#'
#' @description Reads a MIMARKS-compliant (or custom) metadata file in TSV
#'   or CSV format. Optionally matches metadata rows to sample IDs found in
#'   a sequence table and reports matching statistics. Auto-detects standard
#'   MIMARKS fields and logs which are present.
#'
#' @param path Character. Path to a metadata file (TSV or CSV). The format
#'   is auto-detected from the file extension; files with \code{.csv}
#'   extension are read as comma-separated, all others as tab-separated.
#' @param seqtab An optional \code{data.table} in long format with a
#'   \code{sample} column. When provided, metadata rows are matched to
#'   sequence-table sample IDs and only matched rows are returned.
#' @param idColumn Character. Name of the column in the metadata file that
#'   contains sample identifiers (default \code{"sample_name"}).
#'
#' @return A \code{data.table} of metadata. When \code{seqtab} is provided,
#'   only rows matching sequence-table sample IDs are returned and
#'   a \code{matched} attribute records matching statistics.
#'
#' @export
#'
#' @examples
#' ## Create a temporary metadata file
#' tmp <- tempfile(fileext = ".tsv")
#' writeLines(c("sample_name\tsite\tdepth",
#'              "S1\tA\t10",
#'              "S2\tB\t20"), tmp)
#' meta <- loadMetadata(tmp)
#' meta
loadMetadata <- function(path, seqtab = NULL, idColumn = "sample_name") {

    if (!file.exists(path)) {
        stop("Metadata file not found: ", path)
    }

    ## Auto-detect format by extension
    ext <- tolower(tools::file_ext(path))

    if (ext == "csv") {
        message("[loadMetadata] Reading metadata as CSV: ", path)
        meta <- data.table::fread(path, sep = ",", header = TRUE)
    } else if (ext %in% c("tsv", "txt", "tab")) {
        message("[loadMetadata] Reading metadata as TSV: ", path)
        meta <- data.table::fread(path, sep = "\t", header = TRUE)
    } else {
        message("[loadMetadata] Unknown extension '", ext,
                "', trying TSV first")
        meta <- tryCatch(
            data.table::fread(path, sep = "\t", header = TRUE),
            error = function(e) {
                message("[loadMetadata] TSV parse failed, trying CSV")
                data.table::fread(path, sep = ",", header = TRUE)
            }
        )
    }

    message("[loadMetadata] Metadata has ", nrow(meta), " rows and ",
            ncol(meta), " columns")

    ## Check for MIMARKS fields
    mimarks_fields <- c("sample_name", "collection_date", "geo_loc_name",
                         "lat_lon", "depth", "env_broad_scale",
                         "env_local_scale", "env_medium")

    present <- mimarks_fields[mimarks_fields %in% colnames(meta)]
    absent  <- mimarks_fields[!mimarks_fields %in% colnames(meta)]

    if (length(present) > 0) {
        message("[loadMetadata] MIMARKS fields found: ",
                paste(present, collapse = ", "))
    }
    if (length(absent) > 0) {
        message("[loadMetadata] MIMARKS fields absent: ",
                paste(absent, collapse = ", "))
    }

    ## Match to sequence table if provided
    if (!is.null(seqtab)) {
        if (!data.table::is.data.table(seqtab)) {
            stop("'seqtab' must be a data.table with a 'sample' column")
        }

        if (!idColumn %in% colnames(meta)) {
            stop("Sample ID column '", idColumn,
                 "' not found in metadata. Available columns: ",
                 paste(colnames(meta), collapse = ", "))
        }

        seqtab_samples <- unique(seqtab$sample)
        meta_ids       <- as.character(meta[[idColumn]])

        matched      <- intersect(seqtab_samples, meta_ids)
        in_seq_only  <- setdiff(seqtab_samples, meta_ids)
        in_meta_only <- setdiff(meta_ids, seqtab_samples)

        message("[loadMetadata] Matched: ", length(matched), " samples")
        message("[loadMetadata] In seqtab but not metadata: ",
                length(in_seq_only), " samples")
        message("[loadMetadata] In metadata but not seqtab: ",
                length(in_meta_only), " samples")

        if (length(matched) == 0) {
            stop("No samples matched between metadata and seqtab. ",
                 "Check that the '", idColumn,
                 "' column contains the correct sample IDs.")
        }

        ## Subset to matched rows
        meta <- meta[meta[[idColumn]] %in% matched]

        attr(meta, "matched") <- list(
            matched      = matched,
            in_seq_only  = in_seq_only,
            in_meta_only = in_meta_only
        )
    }

    meta
}
