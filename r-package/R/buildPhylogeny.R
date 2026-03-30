#' @title Build a Phylogenetic Tree from ASV Sequences
#'
#' @description Aligns ASV sequences using \code{DECIPHER::AlignSeqs},
#'   trims ragged alignment ends based on a coverage threshold, computes
#'   a combined substitution + indel distance matrix, and builds a
#'   neighbor-joining tree with \code{ape::nj}.
#'
#' @param sequences Character vector of DNA sequences. Names, if present,
#'   are used as tree tip labels; otherwise labels are generated as
#'   \code{ASV_1}, \code{ASV_2}, etc.
#' @param cpus Integer. Number of processors for
#'   \code{DECIPHER::AlignSeqs} (default 1).
#' @param minCoverage Numeric between 0 and 1. Alignment columns where
#'   fewer than this fraction of sequences have a base (i.e., are not
#'   gaps) are trimmed before distance calculation (default 0.6).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{tree}{An object of class \code{ape::phylo}.}
#'   \item{distMatrix}{A \code{dist} object of combined substitution and
#'     indel distances.}
#'   \item{alignment}{A \code{Biostrings::DNAStringSet} of the
#'     (untrimmed) multiple sequence alignment.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' seqs <- c(
#'     "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
#'     "ATGCATGCATGCATGCATGCATGCTTGCATGCATGCATGCATGCATGCATGC",
#'     "ATGCATGCATGCATGCATGCATGCATGCAAGCATGCATGCATGCATGCATGC",
#'     "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCTTGCATGCATGC"
#' )
#' result <- buildPhylogeny(seqs, cpus = 1L)
#' plot(result$tree)
#' }
buildPhylogeny <- function(sequences, cpus = 1L, minCoverage = 0.6) {

    if (!is.character(sequences) || length(sequences) < 3) {
        stop("'sequences' must be a character vector with at least 3 entries")
    }

    ## Assign tip labels
    if (is.null(names(sequences))) {
        seq_ids <- paste0("ASV_", seq_along(sequences))
        names(sequences) <- seq_ids
    } else {
        seq_ids <- names(sequences)
    }

    message("[buildPhylogeny] Aligning ", length(sequences), " sequences...")

    dna <- Biostrings::DNAStringSet(sequences)
    names(dna) <- seq_ids

    alignment <- DECIPHER::AlignSeqs(dna, processors = cpus, verbose = FALSE)

    message("[buildPhylogeny] Alignment: ", Biostrings::width(alignment)[1],
            " positions")

    ## -----------------------------------------------------------------
    ## Trim low-coverage columns
    ## -----------------------------------------------------------------
    aln_matrix <- as.matrix(alignment)
    coverage   <- apply(aln_matrix, 2, function(x) mean(x != "-"))
    keep_cols  <- coverage >= minCoverage

    if (sum(keep_cols) < 50) {
        message("[buildPhylogeny] Very few positions pass coverage filter, ",
                "using full alignment")
        aln_trimmed <- ape::as.DNAbin(alignment)
    } else {
        message("[buildPhylogeny] Trimming alignment: ",
                sum(!keep_cols), " low-coverage positions removed, ",
                sum(keep_cols), " retained")
        trimmed_seqs <- apply(aln_matrix[, keep_cols, drop = FALSE], 1,
                              paste, collapse = "")
        aln_trimmed <- ape::as.DNAbin(
            Biostrings::DNAStringSet(trimmed_seqs)
        )
    }

    ## -----------------------------------------------------------------
    ## Distance matrix and NJ tree
    ## -----------------------------------------------------------------
    message("[buildPhylogeny] Computing distance matrix...")
    dist_raw   <- ape::dist.dna(aln_trimmed, model = "raw",
                                pairwise.deletion = TRUE)
    dist_indel <- ape::dist.dna(aln_trimmed, model = "indel")

    dist_combined <- dist_raw + dist_indel

    ## Replace NaN with max observed distance
    dist_mat <- as.matrix(dist_combined)
    nan_mask <- is.nan(dist_mat)
    if (any(nan_mask)) {
        dist_mat[nan_mask] <- max(dist_mat[!nan_mask], na.rm = TRUE)
        dist_combined <- stats::as.dist(dist_mat)
    }

    message("[buildPhylogeny] Building neighbor-joining tree...")
    tree <- ape::nj(dist_combined)

    message("[buildPhylogeny] Tree: ", ape::Ntip(tree), " tips, ",
            ape::Nnode(tree), " internal nodes")

    list(
        tree       = tree,
        distMatrix = dist_combined,
        alignment  = alignment
    )
}
