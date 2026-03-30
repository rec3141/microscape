test_that("filterSeqtab returns correct structure", {
    library(data.table)

    ## Build a small test dataset
    dt <- data.table(
        sample   = rep(paste0("S", 1:5), each = 6),
        sequence = rep(c(
            paste0(strrep("A", 100), "1"),
            paste0(strrep("A", 100), "2"),
            paste0(strrep("A", 100), "3"),
            paste0(strrep("A", 30), "4"),   # short sequence
            paste0(strrep("A", 100), "5"),
            paste0(strrep("A", 100), "6")
        ), 5),
        count = c(
            100, 50, 30, 5, 1, 200,
            120, 60, 25, 3, 0, 180,
            90,  40, 35, 4, 0, 210,
            110, 55, 28, 6, 0, 190,
            105, 45, 32, 2, 0, 195
        )
    )

    result <- filterSeqtab(dt, minLength = 50, minSamples = 2,
                            minSeqs = 2, minReads = 10)

    ## Check that result is a list with expected names
    expect_type(result, "list")
    expect_named(result, c("filtered", "orphans", "smallSamples", "stats"))

    ## Check that filtered output is a data.table
    expect_s3_class(result$filtered, "data.table")
    expect_true(all(c("sample", "sequence", "count") %in%
                        colnames(result$filtered)))

    ## Check that short sequences were removed
    seq_lengths <- nchar(unique(result$filtered$sequence))
    expect_true(all(seq_lengths >= 50))

    ## Check stats structure
    expect_s3_class(result$stats, "data.frame")
    expect_true("step" %in% colnames(result$stats))
    expect_equal(nrow(result$stats), 5)
})

test_that("filterSeqtab rejects non-data.table input", {
    expect_error(filterSeqtab(data.frame(x = 1)),
                 "must be a data.table")
})

test_that("filterSeqtab rejects missing columns", {
    library(data.table)
    dt <- data.table(sample = "S1", seq = "AAAA", count = 10)
    expect_error(filterSeqtab(dt), "Missing required columns")
})
