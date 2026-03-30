#' @title SparCC-Style Correlation Network
#'
#' @description Computes a compositional correlation network from a
#'   long-format count table. Uses a centered log-ratio (CLR) transform
#'   followed by Pearson correlation to approximate SparCC-style
#'   correlation estimation. Low-prevalence ASVs are filtered out, and
#'   the result is returned as a data.table edge list.
#'
#' @param dt A \code{data.table} in long format with columns \code{sample},
#'   \code{sequence}, and \code{count}.
#' @param minPrevalence Numeric between 0 and 1. Minimum fraction of
#'   samples in which an ASV must be present (with count > 0) to be
#'   included in the analysis (default 0.1).
#' @param minCorrelation Numeric. Minimum absolute correlation for an
#'   edge to be retained in the output (default 0.1).
#'
#' @return A \code{data.table} edge list with columns:
#' \describe{
#'   \item{node1}{Character. First ASV in the pair.}
#'   \item{node2}{Character. Second ASV in the pair.}
#'   \item{correlation}{Numeric. CLR-Pearson correlation.}
#'   \item{weight}{Numeric. Normalized \code{abs(correlation)^3}.}
#'   \item{color}{Character. \code{"blue"} for positive, \code{"red"}
#'     for negative correlations.}
#' }
#'
#' @details
#' The centered log-ratio (CLR) transform addresses the compositionality
#' of sequencing count data. A small pseudocount (0.5) is added before
#' log transformation to handle zeros. This approach avoids the spurious
#' negative correlations that arise when standard Pearson or Spearman
#' correlations are applied to relative abundances.
#'
#' @export
#'
#' @examples
#' library(data.table)
#' set.seed(42)
#' dt <- data.table(
#'     sample   = rep(paste0("S", 1:20), each = 10),
#'     sequence = rep(paste0("seq", 1:10), 20),
#'     count    = rpois(200, lambda = 50)
#' )
#' edges <- sparccNetwork(dt, minPrevalence = 0.1, minCorrelation = 0.1)
#' head(edges)
sparccNetwork <- function(dt, minPrevalence = 0.1, minCorrelation = 0.1) {

    ## Validate input
    if (!data.table::is.data.table(dt)) {
        stop("'dt' must be a data.table with columns: sample, sequence, count")
    }
    required <- c("sample", "sequence", "count")
    missing_cols <- setdiff(required, colnames(dt))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    n_total_asvs <- data.table::uniqueN(dt$sequence)
    n_samples    <- data.table::uniqueN(dt$sample)

    message("[sparccNetwork] Input: ", n_samples, " samples, ",
            n_total_asvs, " ASVs")

    ## -----------------------------------------------------------------
    ## Filter by prevalence
    ## -----------------------------------------------------------------
    min_prev_count <- ceiling(minPrevalence * n_samples)
    asv_prev  <- dt[count > 0, .(prevalence = data.table::uniqueN(sample)),
                    by = sequence]
    keep_seqs <- asv_prev[prevalence >= min_prev_count, sequence]
    n_kept    <- length(keep_seqs)

    message("[sparccNetwork] Prevalence filter (>= ", min_prev_count,
            " samples): ", n_kept, " of ", n_total_asvs, " ASVs retained")

    if (n_kept < 3) {
        stop("Too few ASVs (", n_kept,
             ") pass the prevalence filter. Lower minPrevalence.")
    }

    dt_sub <- dt[sequence %in% keep_seqs]

    ## -----------------------------------------------------------------
    ## Cast to wide matrix
    ## -----------------------------------------------------------------
    dt_wide    <- data.table::dcast(dt_sub, sample ~ sequence,
                                    value.var = "count",
                                    fill = 0L, fun.aggregate = sum)
    sample_ids <- dt_wide$sample
    dt_wide[, sample := NULL]
    count_mat <- as.matrix(dt_wide)
    rownames(count_mat) <- sample_ids

    message("[sparccNetwork] Count matrix: ", nrow(count_mat), " samples x ",
            ncol(count_mat), " ASVs")

    ## -----------------------------------------------------------------
    ## CLR transform + Pearson correlation
    ##
    ## The centered log-ratio (CLR) is defined as:
    ##   clr(x_i) = log(x_i / geometric_mean(x))
    ## A pseudocount of 0.5 handles zeros.
    ## -----------------------------------------------------------------
    pseudo_mat <- count_mat + 0.5
    log_mat    <- log(pseudo_mat)
    geo_means  <- rowMeans(log_mat)
    clr_mat    <- log_mat - geo_means

    cor_mat <- stats::cor(clr_mat, method = "pearson")

    message("[sparccNetwork] Correlation matrix: ",
            nrow(cor_mat), " x ", ncol(cor_mat))

    ## -----------------------------------------------------------------
    ## Melt to edge list (upper triangle only)
    ## -----------------------------------------------------------------
    cor_dt <- data.table::as.data.table(cor_mat, keep.rownames = "node1")
    cor_long <- data.table::melt(cor_dt, id.vars = "node1",
                                  variable.name = "node2",
                                  value.name = "correlation")
    cor_long[, node2 := as.character(node2)]

    ## Upper triangle only
    cor_long <- cor_long[node1 < node2]

    ## Filter by absolute correlation
    cor_long <- cor_long[abs(correlation) > minCorrelation]

    message("[sparccNetwork] Edges with |correlation| > ", minCorrelation,
            ": ", nrow(cor_long))

    ## -----------------------------------------------------------------
    ## Compute edge weights and colors
    ## -----------------------------------------------------------------
    cor_long[, weight := abs(correlation)^3]

    max_weight <- max(cor_long$weight, na.rm = TRUE)
    if (max_weight > 0) {
        cor_long[, weight := weight / max_weight]
    }

    cor_long[, color := ifelse(correlation > 0, "blue", "red")]

    ## Summary
    n_edges    <- nrow(cor_long)
    n_positive <- sum(cor_long$correlation > 0)
    n_negative <- sum(cor_long$correlation < 0)

    message("[sparccNetwork] Network: ", n_edges, " edges (",
            n_positive, " positive, ", n_negative, " negative)")

    cor_long
}
