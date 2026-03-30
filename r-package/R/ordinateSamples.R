#' @title Ordinate Samples and ASVs by Bray-Curtis Distance
#'
#' @description Computes a Bray-Curtis (or other) distance matrix from a
#'   long-format sequence table, then reduces dimensionality via PCA
#'   followed by t-SNE (or PCA alone). Ordination is computed for both
#'   samples and ASVs.
#'
#' @param dt A \code{data.table} in long format with columns \code{sample},
#'   \code{sequence}, and \code{count}.
#' @param method Character. Ordination method: \code{"tsne"} (PCA + t-SNE,
#'   default) or \code{"pca"} (PCA only).
#' @param metric Character. Distance metric passed to \code{stats::dist} or
#'   computed manually. Currently supports \code{"bray"} (Bray-Curtis,
#'   default) and \code{"euclidean"}.
#' @param perplexity Numeric. Perplexity parameter for t-SNE. Automatically
#'   capped at \code{(n - 1) / 3} when the dataset is small (default 30).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{sampleCoords}{A \code{data.frame} with columns \code{label},
#'     \code{Axis1}, \code{Axis2} for sample ordination.}
#'   \item{asvCoords}{A \code{data.frame} with columns \code{label},
#'     \code{Axis1}, \code{Axis2} for ASV ordination.}
#'   \item{distances}{A list with elements \code{samples} and \code{asvs},
#'     each a \code{dist} object.}
#' }
#'
#' @export
#'
#' @examples
#' library(data.table)
#' dt <- data.table(
#'     sample   = rep(paste0("S", 1:6), each = 10),
#'     sequence = rep(paste0("seq", 1:10), 6),
#'     count    = rpois(60, lambda = 50)
#' )
#' res <- ordinateSamples(dt, method = "pca", metric = "bray")
#' head(res$sampleCoords)
ordinateSamples <- function(dt, method = "tsne", metric = "bray",
                             perplexity = 30) {

    ## Validate input
    if (!data.table::is.data.table(dt)) {
        stop("'dt' must be a data.table with columns: sample, sequence, count")
    }
    required <- c("sample", "sequence", "count")
    missing_cols <- setdiff(required, colnames(dt))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    method <- match.arg(method, c("tsne", "pca"))
    metric <- match.arg(metric, c("bray", "euclidean"))

    message("[ordinateSamples] Input: ",
            data.table::uniqueN(dt$sample), " samples, ",
            data.table::uniqueN(dt$sequence), " ASVs")

    ## -----------------------------------------------------------------
    ## Cast to wide proportional matrix
    ## -----------------------------------------------------------------
    dt_wide    <- data.table::dcast(dt, sample ~ sequence,
                                    value.var = "count",
                                    fill = 0L, fun.aggregate = sum)
    sample_ids <- dt_wide$sample
    dt_wide[, sample := NULL]
    mat <- as.matrix(dt_wide)
    rownames(mat) <- sample_ids

    ## Row-wise proportions
    row_totals <- rowSums(mat)
    row_totals[row_totals == 0] <- 1
    prop_mat <- mat / row_totals

    ## -----------------------------------------------------------------
    ## Bray-Curtis distance (manual implementation to avoid external deps)
    ## -----------------------------------------------------------------
    bray_curtis <- function(x) {
        n <- nrow(x)
        d <- matrix(0, n, n)
        for (i in seq_len(n - 1)) {
            for (j in (i + 1):n) {
                num   <- sum(abs(x[i, ] - x[j, ]))
                denom <- sum(x[i, ] + x[j, ])
                d[i, j] <- if (denom > 0) num / denom else 0
                d[j, i] <- d[i, j]
            }
        }
        rownames(d) <- colnames(d) <- rownames(x)
        stats::as.dist(d)
    }

    compute_dist <- function(x) {
        if (metric == "bray") {
            bray_curtis(x)
        } else {
            stats::dist(x, method = "euclidean")
        }
    }

    ## -----------------------------------------------------------------
    ## Ordination helper
    ## -----------------------------------------------------------------
    ordinate <- function(dist_obj, labels, what = "items") {
        dist_mat <- as.matrix(dist_obj)

        if (method == "pca") {
            n_pcs <- min(2, nrow(dist_mat) - 1)
            pca_res <- stats::prcomp(dist_mat, center = TRUE, scale. = FALSE)
            coords <- pca_res$x[, seq_len(n_pcs), drop = FALSE]
            data.frame(
                label = labels,
                Axis1 = coords[, 1],
                Axis2 = if (ncol(coords) >= 2) coords[, 2] else 0,
                stringsAsFactors = FALSE
            )
        } else {
            ## PCA first, then t-SNE
            n_pcs <- min(50, nrow(dist_mat) - 1)
            pca_res <- stats::prcomp(dist_mat, center = TRUE, scale. = FALSE)
            pca_coords <- pca_res$x[, seq_len(n_pcs), drop = FALSE]

            ## Deduplicate
            unique_coords <- unique(pca_coords)
            perp <- min(perplexity, floor((nrow(unique_coords) - 1) / 3))
            if (perp < 1) perp <- 1

            message("[ordinateSamples] Running t-SNE for ", what,
                    " (perplexity = ", perp, ")...")

            ## Use simple PCA-based embedding if Rtsne is not available
            if (requireNamespace("Rtsne", quietly = TRUE)) {
                tsne_res <- Rtsne::Rtsne(unique_coords, dims = 2,
                                          perplexity = perp, theta = 0.5,
                                          pca = FALSE, verbose = FALSE,
                                          check_duplicates = FALSE)
                ## Map back
                unique_key <- apply(unique_coords, 1, paste, collapse = "_")
                orig_key   <- apply(pca_coords, 1, paste, collapse = "_")
                map_idx    <- match(orig_key, unique_key)

                data.frame(
                    label = labels,
                    Axis1 = tsne_res$Y[map_idx, 1],
                    Axis2 = tsne_res$Y[map_idx, 2],
                    stringsAsFactors = FALSE
                )
            } else {
                message("[ordinateSamples] Rtsne not available, ",
                        "falling back to PCA")
                data.frame(
                    label = labels,
                    Axis1 = pca_coords[, 1],
                    Axis2 = if (ncol(pca_coords) >= 2) pca_coords[, 2] else 0,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    ## -----------------------------------------------------------------
    ## Sample ordination
    ## -----------------------------------------------------------------
    message("[ordinateSamples] Computing distances for samples...")
    sample_dist  <- compute_dist(prop_mat)
    sample_coords <- ordinate(sample_dist, labels = rownames(prop_mat),
                              what = "samples")

    ## -----------------------------------------------------------------
    ## ASV ordination (transpose)
    ## -----------------------------------------------------------------
    message("[ordinateSamples] Computing distances for ASVs...")
    asv_mat  <- t(prop_mat)
    asv_dist <- compute_dist(asv_mat)
    asv_coords <- ordinate(asv_dist, labels = rownames(asv_mat),
                           what = "ASVs")

    message("[ordinateSamples] Ordination complete")

    list(
        sampleCoords = sample_coords,
        asvCoords    = asv_coords,
        distances    = list(samples = sample_dist, asvs = asv_dist)
    )
}
