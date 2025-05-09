
test_that("PICBbuild throws error when IN.ALIGNMENTS is NULL", {
    expect_error(
        PICBbuild(
            IN.ALIGNMENTS = NULL,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 0
        ),
        "Please provide IN.ALIGNMENTS !"
    )
})

test_that("PICBbuild throws error when REFERENCE.GENOME is NULL", {
    expect_error(
        PICBbuild(
            IN.ALIGNMENTS = test_alignments,
            REFERENCE.GENOME = NULL,
            VERBOSITY = 0
        ), "Please provide REFERENCE.GENOME !"
    ) 
})


test_that("PICBbuild throws error when IN.ALIGNMENTS is missing 'unique' column", {
    test_alignments_unique_missing <- test_alignments
    test_alignments_unique_missing$unique <- NULL
    expect_error(
        PICBbuild(
            IN.ALIGNMENTS = test_alignments_unique_missing,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6"
        ),
        regexp = "^IN.ALIGNMENTS must contain 'unique' column! "
    )
})

test_that("PICBbuild throws error when IN.ALIGNMENTS is missing 'multi.primary' column", {
    test_alignments_multi_primary_missing <- test_alignments
    test_alignments_multi_primary_missing$multi.primary <- NULL
    expect_error(
        PICBbuild(
            IN.ALIGNMENTS = test_alignments_multi_primary_missing,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 0
        ),
        regexp = "^IN.ALIGNMENTS must contain 'multi.primary' column! "
    )
})

test_that("PICBbuild warns when IN.ALIGNMENTS is missing 'multi.secondary' column", {
    test_alignments_multi_secondary_missing <- test_alignments
    test_alignments_multi_secondary_missing$multi.secondary <- NULL
    expect_warning(
        PICBbuild(
            IN.ALIGNMENTS = test_alignments_multi_secondary_missing,
            REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
            VERBOSITY = 0
        ),
        regexp = "^IN.ALIGNMENTS does not contain secondary multimappers"
    )
})

test_that("PICBbuild builds piRNA seeds, cores, and clusters with default parameters correctly", {
    expect_type(test_ranges, "list")
    expect_named(test_ranges, c("seeds", "cores", "clusters"))
    expect_true(all(sapply(test_ranges, function(x) inherits(x, "GRanges"))))

    expect_equal(length(test_ranges$seeds), 42)
    expect_equal(length(test_ranges$cores), 38)
    expect_equal(length(test_ranges$clusters), 27)
    # Check metadata columns
    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
    expect_equal(colnames(mcols(test_ranges$seeds)), expectedOutputCols)
    expect_equal(colnames(mcols(test_ranges$cores)), expectedOutputCols)
    expect_equal(colnames(mcols(test_ranges$clusters)), c(expectedOutputCols, "type"))
})


test_that("PICBbuild example subset builds piRNA seeds, cores, and clusters correctly with defined LIBRARY.SIZE", {
    result <- PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        LIBRARY.SIZE = 12799826
    )
    expect_type(result, "list")
    expect_named(result, c("seeds", "cores", "clusters"))
    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))

    expect_equal(length(result$seeds), 25)
    expect_equal(length(result$cores), 18)
    expect_equal(length(result$clusters), 4)
    # Check metadata columns
    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
    expect_equal(colnames(mcols(result$seeds)), expectedOutputCols)
    expect_equal(colnames(mcols(result$cores)), expectedOutputCols)
    expect_equal(colnames(mcols(result$clusters)), c(expectedOutputCols, "type"))
})

#optional additional testing
#test_that("PICBbuild handles custom sliding window sizes correctly", {
#    custom_uniq_window_width <- 500
#    custom_uniq_window_step <- 50
#    custom_primary_mult_window_width <- 400
#    custom_primary_mult_window_step <- 40
#    custom_secondary_mult_window_width <- 1200
#    custom_secondary_mult_window_step <- 120
#    result <- PICBbuild(
#        IN.ALIGNMENTS = test_alignments,
#        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#        UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH = custom_uniq_window_width,
#        UNIQUEMAPPERS.SLIDING.WINDOW.STEP = custom_uniq_window_step,
#        PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH = custom_primary_mult_window_width,
#        PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.STEP = custom_primary_mult_window_step,
#        SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH = custom_secondary_mult_window_width,
#        SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.STEP = custom_secondary_mult_window_step
#    )
#    expect_type(result, "list")
#    expect_named(result, c("seeds", "cores", "clusters"))
#    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))
#
#    expect_equal(length(result$seeds), 43)
#    expect_equal(length(result$cores), 36)
#    expect_equal(length(result$clusters), 28)
#    # Check metadata columns
#    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
#    expect_equal(colnames(mcols(result$seeds)), expectedOutputCols)
#    expect_equal(colnames(mcols(result$cores)), expectedOutputCols)
#    expect_equal(colnames(mcols(result$clusters)), c(expectedOutputCols, "type"))
#})


#test_that("PICBbuild filters based on custom thresholds correctly", {
#    custom_min_overlap <- 10
#    custom_threshold_seeds_gap <- 100
#    custom_threshold_cores_gap <- 200
#    custom_threshold_clusters_gap <- 300
#    result <- PICBbuild(
#        IN.ALIGNMENTS = test_alignments,
#        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#        MIN.OVERLAP = custom_min_overlap,
#        THRESHOLD.SEEDS.GAP = custom_threshold_seeds_gap,
#        THRESHOLD.CORES.GAP = custom_threshold_cores_gap,
#        THRESHOLD.CLUSTERS.GAP = custom_threshold_clusters_gap
#    )
#    expect_type(result, "list")
#    expect_named(result, c("seeds", "cores", "clusters"))
#    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))
#
#
#    expect_equal(length(result$seeds), 41)
#    expect_equal(length(result$cores), 31)
#    expect_equal(length(result$clusters), 20)
#    # check metacolumns
#    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
#    expect_equal(colnames(mcols(result$seeds)), expectedOutputCols)
#    expect_equal(colnames(mcols(result$cores)), expectedOutputCols)
#    expect_equal(colnames(mcols(result$clusters)), c(expectedOutputCols, "type"))
#})


#test_that("PICBbuild provides non-normalized statistics when PROVIDE.NON_NORMALIZED = TRUE", {
#    result <- PICBbuild(
#        IN.ALIGNMENTS = test_alignments,
#        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#        PROVIDE.NON.NORMALIZED = TRUE
#    )
#    expect_type(result, "list")
#    expect_named(result, c("seeds", "cores", "clusters"))
#
#   # Check metadata columns
#    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
#    addNonNormOutputCols <- c("uniq_reads", "multimapping_reads_primary_alignments", "all_reads_primary_alignments", "width_covered_by_unique_alignments", "uniq_sequences")
#    expect_equal(length(setdiff(colnames(mcols(result$seeds)), c(expectedOutputCols, addNonNormOutputCols))), 0)
#    expect_equal(length(setdiff(colnames(mcols(result$cores)), c(expectedOutputCols, addNonNormOutputCols))), 0)
#    expect_equal(length(setdiff(colnames(mcols(result$seeds)), c(expectedOutputCols, addNonNormOutputCols, "type"))), 0)
#})

test_that("PICBbuild loads BAM file with default parameters the same with using the same reference genome but different methods", {
    # Using BSgenome as myGenome shown above. 

    # Using chromosome names and lengths
    seqnames <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
    seqlengths <- c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352)
    myGenome <- Seqinfo(seqnames = seqnames, seqlengths = seqlengths)

    result <- PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = myGenome
    )
    expect_equal(result, test_ranges)

    # Using genome name
    myGenome <- GenomeInfoDb::Seqinfo(genome = "dm6")
    result <- PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = myGenome,
        VERBOSITY = 0
    )
    expect_equal(result, test_ranges)

    # Using a FASTA file
    chr2L_seq <- BSgenome.Dmelanogaster.UCSC.dm6[["chr2L"]]
    chr2L_seq_set <- DNAStringSet(chr2L_seq)
    names(chr2L_seq_set) <- "chr2L"
    temp_fasta <- tempfile(fileext = ".fasta")
    writeXStringSet(chr2L_seq_set, temp_fasta)
    myGenome <- PICBloadfasta(FASTA.NAME = temp_fasta)
    result <- PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = myGenome,
        VERBOSITY = 0
    )
    # seqinfo is expected to be different (all standard chromosomes in above examples vs just chr2L in subset genome)
    seqInfoExp <- seqinfo(test_ranges$seeds)
    seqinfo(result$seeds) <- seqInfoExp
    seqinfo(result$cores) <- seqInfoExp
    seqinfo(result$clusters) <- seqInfoExp

    expect_equal(result, test_ranges)
    unlink(temp_fasta)
})


test_that("PICBbuild changes seqlevels style when SEQ.LEVELS.STYLE is specified", {
    custom_style <- "NCBI"
    result <- PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        SEQ.LEVELS.STYLE = custom_style,
        VERBOSITY = 0
    )
    expect_type(result, "list")
    expect_named(result, c("seeds", "cores", "clusters"))
    # Check if seqlevels style is updated
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$seeds)[1], custom_style)
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$cores)[1], custom_style)
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$clusters)[1], custom_style)
})

test_that("PICBbuild computes 1U and 10A fractions when COMPUTE.1U.10A.FRACTIONS = TRUE", {
    # for error when column seq not provided in PICBload() output (test_alignments)
    expect_error(
        PICBbuild(
        IN.ALIGNMENTS = test_alignments,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        COMPUTE.1U.10A.FRACTIONS = TRUE
        ),
        regexp = "^Alignments 'unique' does not contain sequence"
    )
    # proper calculations with column seq in test_alignments
    ouputPICBloadWseq <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        GET.ORIGINAL.SEQUENCE = TRUE,
        VERBOSE = FALSE
    )

    result <- PICBbuild(
        IN.ALIGNMENTS = ouputPICBloadWseq,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        COMPUTE.1U.10A.FRACTIONS = TRUE
    )

    expect_type(result, "list")
    expect_named(result, c("seeds", "cores", "clusters"))

    # Check that the fraction columns are present in metadata
    expectedOutputCols <- c("width_in_nt", "uniq_reads_FPM", "multimapping_reads_primary_alignments_FPM", "all_reads_primary_alignments_FPM", "uniq_reads_FPKM", "multimapping_reads_primary_alignments_FPKM", "all_reads_primary_alignments_FPKM", "fraction_of_width_covered_by_unique_alignments")
    fraction1U10ACols <- c("oneU.frac.unique", "tenA.frac.unique", "oneU.frac.multi.primary", "tenA.frac.multi.primary", "oneU.frac.multi.secondary", "tenA.frac.multi.secondary")
    expect_equal(length(setdiff(colnames(mcols(result$seeds)), c(expectedOutputCols, fraction1U10ACols))), 0)
    expect_equal(length(setdiff(colnames(mcols(result$cores)), c(expectedOutputCols, fraction1U10ACols))), 0)
    expect_equal(length(setdiff(colnames(mcols(result$seeds)), c(expectedOutputCols, fraction1U10ACols, "type"))), 0)
})

test_that("PICBbuild handles empty alignments gracefully", {
    # Create empty GRanges objects
    empty_gr <- GRanges()
    alignments_empty <- list(
        unique = empty_gr,
        multi.primary = empty_gr,
        multi.secondary = empty_gr
    )
    expect_error(
        PICBbuild(
        IN.ALIGNMENTS = alignments_empty,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6"
        ), "The IN.ALIGNMENTS must contain NH information !"
    )
})

test_that("PICBbuild VERBOSITY parameter works correctly", {
  
  # Test VERBOSITY = 0 (silent)
  silent_messages <- capture_messages({
    PICBbuild(
      IN.ALIGNMENTS = test_alignments,
      REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
      VERBOSITY = 0
    )
  })
  
  # Test VERBOSITY = 1 (basic messages)
  basic_messages <- capture_messages({
    PICBbuild(
      IN.ALIGNMENTS = test_alignments,
      REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
      VERBOSITY = 1
    )
  })
  
  # Test VERBOSITY = 2 (detailed messages)
  detailed_messages <- capture_messages({
    PICBbuild(
      IN.ALIGNMENTS = test_alignments,
      REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
      VERBOSITY = 3
    )
  })
  
  # Test expectations
  expect_identical(length(silent_messages), 0L, "VERBOSITY=0 should not produce any messages")
  
  expect_true(length(basic_messages) == 4, "VERBOSITY=1 should produce messages")
  expect_true(
    any(grepl("Starting", basic_messages)) && 
    any(grepl("Done", basic_messages)),
    "VERBOSITY=1 should contain basic progress messages"
  )
  
  expect_true(length(detailed_messages) > length(basic_messages), 
              "VERBOSITY=2 should produce more messages than VERBOSITY=1")
  expect_true(
    any(grepl("Sliding window analysis", detailed_messages)) && 
    any(grepl("Annotating", detailed_messages)),
    "VERBOSITY=2 should contain detailed progress messages"
  )
})