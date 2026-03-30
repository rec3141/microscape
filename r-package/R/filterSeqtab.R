#' @title Filter a Long-Format Sequence Table
#'
#' @description Applies a cascade of quality-control filters to a long-format
#'   sequence table (data.table). Filters are applied in order: (1) minimum
#'   sequence length, (2) minimum sample prevalence, (3) minimum total
#'   abundance, and (4) minimum per-sample read depth. Operates entirely in
#'   long format so memory use is proportional to non-zero entries.
#'
#' @param dt A \code{data.table} in long format with columns \code{sample},
#'   \code{sequence}, and \code{count}.
#' @param minLength Integer. Minimum sequence length in base pairs. ASVs
#'   shorter than this are removed (default 50).
#' @param minSamples Integer. Minimum number of samples an ASV must appear
#'   in. ASVs below this threshold are classified as orphans (default 2).
#' @param minSeqs Integer. Minimum total read count across all samples for
#'   an ASV to be retained (default 2).
#' @param minReads Integer. Minimum total read count for a sample to be
#'   retained (default 1000).
#'
#' @return A list with four elements:
#' \describe{
#'   \item{filtered}{A \code{data.table} containing the filtered sequence
#'     table in long format.}
#'   \item{orphans}{A \code{data.table} of ASVs removed by the prevalence
#'     filter.}
#'   \item{smallSamples}{A \code{data.table} of samples removed by the
#'     depth filter.}
#'   \item{stats}{A \code{data.frame} summarizing ASVs and samples removed
#'     at each filtering step.}
#' }
#'
#' @export
#'
#' @examples
#' library(data.table)
#' dt <- data.table(
#'     sample = rep(paste0("S", 1:5), each = 10),
#'     sequence = rep(paste0(strrep("A", 100), seq_len(10)), 5),
#'     count = sample(1:100, 50, replace = TRUE)
#' )
#' result <- filterSeqtab(dt, minLength = 50, minSamples = 2,
#'                         minSeqs = 2, minReads = 10)
#' names(result)
filterSeqtab <- function(dt, minLength = 50L, minSamples = 2L,
                          minSeqs = 2L, minReads = 1000L) {

    ## Validate input
    if (!data.table::is.data.table(dt)) {
        stop("'dt' must be a data.table with columns: sample, sequence, count")
    }
    required <- c("sample", "sequence", "count")
    missing_cols <- setdiff(required, colnames(dt))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    ## Work on a copy to avoid modifying the input
    dt <- data.table::copy(dt)

    total_input  <- sum(dt$count)
    n_input_asvs <- data.table::uniqueN(dt$sequence)

    message("[filterSeqtab] Input: ", data.table::uniqueN(dt$sample),
            " samples, ", n_input_asvs, " ASVs, ", total_input, " reads")

    ## 1. Remove short sequences
    seq_lengths <- dt[, .(len = nchar(sequence[1])), by = sequence]
    short_seqs  <- seq_lengths[len < minLength, sequence]
    n_short     <- length(short_seqs)

    if (n_short > 0) {
        message("[filterSeqtab] Removing ", n_short, " ASVs shorter than ",
                minLength, " bp")
        dt <- dt[!sequence %in% short_seqs]
    }

    ## 2. Remove low-prevalence ASVs (orphans)
    seq_prevalence <- dt[, .(n_samples = data.table::uniqueN(sample)),
                         by = sequence]
    orphan_seqs    <- seq_prevalence[n_samples < minSamples, sequence]
    n_orphans      <- length(orphan_seqs)

    dt_orphans <- dt[sequence %in% orphan_seqs]
    dt         <- dt[!sequence %in% orphan_seqs]

    message("[filterSeqtab] Removed ", n_orphans,
            " orphan ASVs (present in < ", minSamples, " samples)")

    ## 3. Remove low-abundance ASVs
    seq_abundance <- dt[, .(total = sum(count)), by = sequence]
    rare_seqs     <- seq_abundance[total < minSeqs, sequence]
    n_rare        <- length(rare_seqs)

    dt <- dt[!sequence %in% rare_seqs]

    message("[filterSeqtab] Removed ", n_rare, " rare ASVs (< ", minSeqs,
            " total reads)")

    ## 4. Remove shallow samples
    sample_depth  <- dt[, .(total = sum(count)), by = sample]
    small_samples <- sample_depth[total < minReads, sample]
    n_small       <- length(small_samples)

    dt_small <- dt[sample %in% small_samples]
    dt       <- dt[!sample %in% small_samples]

    message("[filterSeqtab] Removed ", n_small, " samples (< ", minReads,
            " reads)")

    ## Summary
    n_final_samples <- data.table::uniqueN(dt$sample)
    n_final_asvs    <- data.table::uniqueN(dt$sequence)
    n_final_reads   <- sum(dt$count)
    pct_retained    <- round(n_final_reads / max(total_input, 1) * 100, 1)

    message("[filterSeqtab] Final: ", n_final_samples, " samples, ",
            n_final_asvs, " ASVs, ", n_final_reads, " reads (",
            pct_retained, "% of input)")

    ## Build stats summary
    stats <- data.frame(
        step               = c("length", "prevalence", "abundance",
                                "depth", "final"),
        asvs_removed       = c(n_short, n_orphans, n_rare, NA, NA),
        samples_removed    = c(NA, NA, NA, n_small, NA),
        remaining_samples  = c(NA, NA, NA, NA, n_final_samples),
        remaining_asvs     = c(NA, NA, NA, NA, n_final_asvs),
        pct_reads_retained = c(NA, NA, NA, NA, pct_retained),
        stringsAsFactors   = FALSE
    )

    list(
        filtered     = dt,
        orphans      = dt_orphans,
        smallSamples = dt_small,
        stats        = stats
    )
}
