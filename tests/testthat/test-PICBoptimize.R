optOutput <- PICBoptimize(
    IN.ALIGNMENTS = test_alignments,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSITY = 0,
    MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2)
)

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


test_that("PICBoptimize returns correct number of rows for MIN.UNIQUE.ALIGNMENTS.PER.WINDOW", {
    expect_equal(nrow(optOutput), 2) 
})

test_that("PICBoptimize handles different verbosity levels", {
    # Test verbosity = 2 (should show version message)
    expect_message(
        PICBoptimize(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 2,
            MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2)
        ),
        regexp = "^PICB v"
    )
})


test_that("PICBoptimize handles basic parameter iteration correctly", {
    # Test single parameter iteration
    result1 <- PICBoptimize(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2, 3)
    )
    expect_equal(nrow(result1), 3)
    expect_true(all(c("MIN.UNIQUE.ALIGNMENTS.PER.WINDOW") %in% colnames(result1)))
    
    # Test multiple parameter iteration
    result2 <- PICBoptimize(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2),
        MIN.SEED.LENGTH = c(800, 1000)
    )
    expect_equal(nrow(result2), 4)  # 2x2 combinations
    expect_true(all(c("MIN.UNIQUE.ALIGNMENTS.PER.WINDOW", "MIN.SEED.LENGTH") %in% colnames(result2)))
})

test_that("PICBoptimize handles library size correctly", {
    custom_library_size <- 1000000
    # Library size warning should be triggered
    expect_warning(
        PICBoptimize(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            LIBRARY.SIZE = custom_library_size,
            MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2)
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
    expect_true(all(expected_cols %in% colnames(optOutput)))
})

test_that("PICBoptimize handles PROVIDE.INFO.SEEDS.AND.CORES correctly", {
    result <- PICBoptimize(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSITY = 0,
        PROVIDE.INFO.SEEDS.AND.CORES = TRUE,
        MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2)
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
        "MIN.UNIQUE.ALIGNMENTS.PER.WINDOW",  # Include any fixed columns first
        outer(base_metrics, region_types, function(x, y) paste0(x, ".", y))
    )
    
    expect_true(all(expected_cols %in% colnames(result)))
})