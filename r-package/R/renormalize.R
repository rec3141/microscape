#' @title Renormalize Counts by Taxonomic Group
#'
#' @description Classifies ASVs into taxonomic groups (prokaryote,
#'   chloroplast, mitochondria, eukaryote, unknown) using a taxonomy matrix,
#'   then normalizes counts to within-group proportions per sample. This
#'   prevents organelle or cross-domain reads from distorting relative
#'   abundance estimates.
#'
#' @param dt A \code{data.table} in long format with columns \code{sample},
#'   \code{sequence}, and \code{count}.
#' @param taxa A character matrix of taxonomy assignments with rows named by
#'   sequence and columns for taxonomic ranks (e.g., Domain, Phylum, ...,
#'   Genus). SILVA-style rank names are expected.
#'
#' @return A named list of \code{data.table}s, one per taxonomic group.
#'   Each table has columns \code{sample}, \code{sequence}, \code{count},
#'   and \code{proportion} (within-group per-sample proportion).
#'
#' @export
#'
#' @examples
#' library(data.table)
#' dt <- data.table(
#'     sample   = rep(c("S1", "S2"), each = 3),
#'     sequence = rep(c("AAAA", "CCCC", "GGGG"), 2),
#'     count    = c(100, 50, 10, 80, 60, 20)
#' )
#' taxa <- matrix(c("Bacteria", "Cyanobacteria", "Cyanobacteriia",
#'                   "Chloroplast", NA, NA,
#'                   "Bacteria", "Proteobacteria", "Alphaproteobacteria",
#'                   "Rickettsiales", "Mitochondria", NA,
#'                   "Bacteria", "Firmicutes", "Bacilli",
#'                   "Lactobacillales", "Lactobacillaceae", "Lactobacillus"),
#'                nrow = 3, byrow = TRUE,
#'                dimnames = list(c("AAAA", "CCCC", "GGGG"),
#'                                c("Domain", "Phylum", "Class",
#'                                  "Order", "Family", "Genus")))
#' groups <- renormalize(dt, taxa)
#' names(groups)
renormalize <- function(dt, taxa) {

    ## Validate inputs
    if (!data.table::is.data.table(dt)) {
        stop("'dt' must be a data.table with columns: sample, sequence, count")
    }
    required <- c("sample", "sequence", "count")
    missing_cols <- setdiff(required, colnames(dt))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    if (!is.matrix(taxa)) {
        stop("'taxa' must be a character matrix with rows named by sequence")
    }

    dt <- data.table::copy(dt)

    message("[renormalize] Input: ", data.table::uniqueN(dt$sample),
            " samples, ", data.table::uniqueN(dt$sequence), " ASVs, ",
            sum(dt$count), " reads")

    ## -----------------------------------------------------------------
    ## Classify ASVs into taxonomic groups
    ## -----------------------------------------------------------------
    all_seqs <- unique(dt$sequence)

    ## Helper: find sequences matching a value at a given rank
    find_at_rank <- function(rank_name, values) {
        if (!rank_name %in% colnames(taxa)) return(character(0))
        rownames(taxa)[taxa[, rank_name] %in% values &
                           !is.na(taxa[, rank_name])]
    }

    chloroplast_seqs  <- find_at_rank("Order", "Chloroplast")
    mitochondria_seqs <- find_at_rank("Family", "Mitochondria")
    eukaryote_seqs    <- find_at_rank("Domain", "Eukaryota")
    prokaryote_seqs   <- find_at_rank("Domain", c("Bacteria", "Archaea"))

    ## Remove organelles from prokaryote set
    prokaryote_seqs <- setdiff(prokaryote_seqs,
                               c(chloroplast_seqs, mitochondria_seqs))

    ## Build group map
    group_map <- data.table::data.table(
        sequence = rownames(taxa), group = "unknown"
    )
    group_map[sequence %in% chloroplast_seqs,  group := "chloroplast"]
    group_map[sequence %in% mitochondria_seqs, group := "mitochondria"]
    group_map[sequence %in% eukaryote_seqs,    group := "eukaryote"]
    group_map[sequence %in% prokaryote_seqs,   group := "prokaryote"]

    ## Sequences in dt but not in taxonomy
    missing_seqs <- setdiff(all_seqs, group_map$sequence)
    if (length(missing_seqs) > 0) {
        group_map <- rbind(
            group_map,
            data.table::data.table(sequence = missing_seqs, group = "unknown")
        )
    }

    ## Log group sizes
    group_counts <- group_map[, .N, by = group][order(-N)]
    for (i in seq_len(nrow(group_counts))) {
        message("[renormalize]   ", group_counts$group[i], ": ",
                group_counts$N[i], " ASVs")
    }

    ## -----------------------------------------------------------------
    ## Merge groups and normalize within each group
    ## -----------------------------------------------------------------
    dt_grouped <- merge(dt, group_map, by = "sequence", all.x = TRUE)
    dt_grouped[is.na(group), group := "unknown"]

    dt_grouped[, total_in_group := sum(count), by = .(sample, group)]
    dt_grouped[, proportion := count /
                   ifelse(total_in_group > 0, total_in_group, 1)]
    dt_grouped[is.nan(proportion), proportion := 0]

    ## -----------------------------------------------------------------
    ## Split into per-group tables
    ## -----------------------------------------------------------------
    groups <- unique(dt_grouped$group)
    group_tables <- list()

    for (g in groups) {
        dt_g <- dt_grouped[group == g,
                           .(sample, sequence, count, proportion)]
        group_tables[[g]] <- dt_g
        message("[renormalize] ", g, ": ",
                data.table::uniqueN(dt_g$sequence), " ASVs, ",
                sum(dt_g$count), " reads across ",
                data.table::uniqueN(dt_g$sample), " samples")
    }

    group_tables
}
