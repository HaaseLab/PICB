test_that("PICBoptimize throws error when IN.ALIGNMENTS is NULL", {
    expect_error(
        PICBoptimize(
            IN.ALIGNMENTS = NULL,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 0
        ),
        "Please provide IN.ALIGNMENTS !"
    )
})

test_that("PICBoptimize throws error when REFERENCE.GENOME is NULL", {
    expect_error(
        PICBoptimize(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = NULL,
            VERBOSITY = 0
        ), "Please provide REFERENCE.GENOME !"
    ) 
})

test_that("PICBoptimize throws error when nothing to optimize for", {
    expect_error(
        PICBoptimize(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 0
        ), "Provide arguments to iterave over. See example."
    ) 
})

test_that("PICBoptimize handles basic parameter iteration correctly", {
    # Test single parameter iteration
    expect_equal(nrow(test_optimization), 2)
    expect_true(all(c("MIN.UNIQUE.ALIGNMENTS.PER.WINDOW") %in% colnames(test_optimization)))
    
    # Test multiple parameter iteration
    result2 <- PICBoptimize(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2),
        THRESHOLD.CLUSTERS.GAP = c(0, 50)
    )
    expect_equal(nrow(result2), 4)  # 2x2 combinations
    expect_true(all(c("MIN.UNIQUE.ALIGNMENTS.PER.WINDOW", "THRESHOLD.CLUSTERS.GAP") %in% colnames(result2)))
})

test_that("PICBoptimize handles library size correctly", {
    custom_library_size <- 1000000
    # Library size warning should be triggered
    expect_warning(
        PICBoptimize(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            LIBRARY.SIZE = custom_library_size,
            THRESHOLD.SEEDS.GAP = c(0, 50)
        ),
        regexp = "^The total number of primary alignments is not equal to the total number of read names"
    )
})

test_that("PICBoptimize output contains expected columns", {
    
    expected_cols <- c(
        "MIN.UNIQUE.ALIGNMENTS.PER.WINDOW",
        "number.of.clusters",
        "total.width.clusters",
        "reads.explained.by.clusters",
        "fraction.of.library.explained.by.clusters",
        "mean.RPKM.clusters",
        "fraction.of.genome.space.clusters"
    )
    expect_true(all(expected_cols %in% colnames(test_optimization)))
})

test_that("PICBoptimize handles PROVIDE.INFO.SEEDS.AND.CORES correctly", {
    result <- PICBoptimize(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        PROVIDE.INFO.SEEDS.AND.CORES = TRUE,
        THRESHOLD.CORES.GAP = c(0, 50)
    )
    
    # Base metrics that apply to each type
    base_metrics <- c(
        "number.of",
        "total.width",
        "reads.explained.by",
        "fraction.of.library.explained.by",
        "mean.RPKM",
        "fraction.of.genome.space"
    )
    
    # Types of regions
    region_types <- c("seeds", "cores", "clusters")
    
    # Generate all expected column combinations
    expected_cols <- c(
        "THRESHOLD.CORES.GAP",  # Include any fixed columns first
        outer(base_metrics, region_types, function(x, y) paste0(x, ".", y))
    )
    
    expect_true(all(expected_cols %in% colnames(result)))
})