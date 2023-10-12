test_that("pbmc_10k() works", {
    x <- pbmc_10x()
    classes <- unlist(lapply(x, class), use.names = FALSE)
    expect_true(all(classes == c(
        "BamFile",
        "GRanges",
        "SingleCellExperiment"
    )))
    expect_equal(length(x$sites), 15638648)
})

test_that("NA12878() works", {
    x <- NA12878()
    fa <- Rsamtools::scanFa(x$fasta)
    expect_equal(length(x$bam), 2)
    expect_equal(length(x$snps), 380175)
    expect_s4_class(fa, "DNAStringSet")
})

test_that("GSE99249() works", {
    x <- GSE99249()
    fa <- Rsamtools::scanFa(x$fasta)
    expect_equal(length(x$bams), 6)
    expect_equal(length(x$sites), 15638648)
    expect_s4_class(fa, "DNAStringSet")
})
