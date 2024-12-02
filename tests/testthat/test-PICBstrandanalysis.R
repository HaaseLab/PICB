# expect_equal(myFunction(testData), expectedOutput)

# invalid IN.ALIGNMENTS
test_that("PICBstrandanalysis stops with invalid IN.ALIGNMENTS", {
    alignments <- list(unique = GenomicRanges::GRanges())
    ranges <- GenomicRanges::GRanges()

    expect_error(PICBstrandanalysis(NULL, ranges), "IN.ALIGNMENTS is either NULL or contains no rows. Please provide valid IN.ALIGNMENTS from PICBload!")
})

# invalid IN.RANGES
test_that("PICBstrandanalysis stops with invalid IN.RANGES", {
    alignments <- list(unique = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1000, 2000), width = 50),
        strand = c("+", "-"),
        NH = c(1, 1)
    ))
    ranges <- GenomicRanges::GRanges()
    expect_error(PICBstrandanalysis(alignments, NULL), "IN.RANGES is either NULL or contains no rows. Use the GRanges object seeds, cores, or clusters from PICBbuild to create valid IN.RANGES!")
    expect_error(PICBstrandanalysis(alignments, ranges), "IN.RANGES is either NULL or contains no rows. Use the GRanges object seeds, cores, or clusters from PICBbuild to create valid IN.RANGES!")
    expect_error(PICBstrandanalysis(alignments, list(ranges)), "IN.RANGES must be a single GRanges object: seeds, cores, or clusters from PICBbuild!")
})

# valid IN.ALIGNMENTS and IN.RANGES expect_equal
test_that("PICBstrandanalysis returns expected output when alignments on both strands", {
    alignments <- list(unique = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(15, 27, 30, 19, 43, 52, 75, 80), width = c(26, 27, 27, 27, 27, 26, 27, 27)),
        strand = c("+", "+", "+", "+", "-", "-", "-", "-"),
        NH = c(1, 1, 1, 1, 1, 1, 1, 1)
    ), multi.primary = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(20, 30, 40, 60, 68, 81), width = c(26, 27, 27, 27, 27, 27)),
        strand = c("+", "-", "-", "-", "+", "-"),
        NH = c(3, 2, 4, 100, 34, 90)
    ))

    ranges <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(10, 40), width = c(40, 60)),
        strand = c("+", "-")
    )
    # piC on plus strand: sense/antisense = 4/1;
    # piC on minus strand: sense/antisense = 4/4;

    # Note: With piRNA input, we anticipate a large number of alignments.
    # To prevent division by zero when calculating the sense-to-antisense ratio (s_as_ratio), we add a pseudocount of 1 to both the sense and antisense counts.
    # Consequently, the expected value in the s_as_ratio column is adjusted to 2.5 and 1.0 instead of 4.0 and 1.0 in this mock dataset.
    # This adjustment ensures numerical stability during computation and does not affect the accuracy of results when working with real piRNA data, where counts are sufficiently large.

    # piC on plus strand: sense+1/antisense+1 = 5/2 = 2.5
    # piC on minus strand: sense+1/antisense+1 = 5/5 = 1.0

    expect_equal(
        PICBstrandanalysis(alignments, ranges),
        GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(10, 40), width = c(40, 60)),
        strand = c("+", "-"),
        s_as_ratio = c(2.5, 1.0)
        )
    )
})


# valid IN.ALIGNMENTS and IN.RANGES but just one strand in IN.RANGES, expect_equal
test_that("PICBstrandanalysis returns expected output when no clusters on plus strand", {
    alignments <- list(unique = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(15, 27, 30, 19, 43, 52, 75, 80), width = c(26, 27, 27, 27, 27, 26, 27, 27)),
        strand = c("+", "+", "+", "+", "-", "-", "-", "-"),
        NH = c(1, 1, 1, 1, 1, 1, 1, 1)
    ), multi.primary = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(20, 30, 40, 60, 68, 81), width = c(26, 27, 27, 27, 27, 27)),
        strand = c("+", "-", "-", "-", "+", "-"),
        NH = c(3, 2, 4, 100, 34, 90)
    ))

    ranges <- GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(10), width = c(40)),
        strand = c("+")
    )

    expect_equal(
        PICBstrandanalysis(alignments, ranges),
        GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(10), width = c(40)),
        strand = c("+"),
        s_as_ratio = c(2.5)
        )
    )

    ranges <- GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(40), width = c(60)),
        strand = c("-")
    )

    expect_equal(
        PICBstrandanalysis(alignments, ranges),
        GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(40), width = c(60)),
        strand = c("-"),
        s_as_ratio = c(1.0)
        )
    )
})

# valid IN.RANGES but IN.ALIGNMENTS not corresponding to IN.RANGES, expect_equal
test_that("PICBstrandanalysis returns expected output when IN.ALIGNMENTS not corresponding to IN.RANGES", {
    alignments <- list(unique = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(15, 27, 30, 19), width = c(26, 27, 27, 27)),
        strand = c("+", "+", "+", "+"),
        NH = c(1, 1, 1, 1)
    ), multi.primary = GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(20, 68), width = c(26, 27)),
        strand = c("+", "+"),
        NH = c(3, 34)
    ))

    ranges <- GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(40), width = c(60)),
        strand = c("-")
    )

    result <- expect_warning(
        PICBstrandanalysis(alignments, ranges),
        "No overlaps detected between IN.ALIGNMENTS and clusters on the minus strand. Please verify that IN.ALIGNMENTS and IN.RANGES are valid and correctly correspond to each other. Continuing ..."
    )
    # piC on minus strand: sense+1/antisense+1 = 1/5 = 0.2
    expect_equal(
        result,
        GenomicRanges::GRanges(
        seqnames = c("chr1"),
        ranges = IRanges::IRanges(start = c(40), width = c(60)),
        strand = c("-"),
        s_as_ratio = c(0.2)
        )
    )
})


# test with example data

test_that("PICBload loads BAM file with default parameters correctly", {
    bam <- system.file("extdata", "Fly_Ov1_chr2L_20To21mb.bam", package = "PICB")
    resultPICBload <- PICBload(
        BAMFILE = bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSE = FALSE
    )
    expect_type(resultPICBload, "list")
    expect_named(resultPICBload, c("unique", "multi.primary", "multi.secondary"))
    expect_true(all(sapply(resultPICBload, function(x) inherits(x, "GRanges"))))

    myClusters <- PICBbuild(
        IN.ALIGNMENTS = resultPICBload,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        LIBRARY.SIZE = 12799826
    )$clusters
    expect_equal(
        length(PICBstrandanalysis(IN.ALIGNMENTS = resultPICBload, IN.RANGES = myClusters)$s_as_ratio), 4
    )
    expect_equal(
        round(PICBstrandanalysis(IN.ALIGNMENTS = resultPICBload, IN.RANGES = myClusters)$s_as_ratio, digits = 3), round(c(3.9041074, 1.6623939, 0.2567568, 0.601542), digits = 3)
    )
})
