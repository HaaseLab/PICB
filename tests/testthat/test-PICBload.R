test_that("example data (bam and corresponding bai file) is present", {
    expect_true(file.exists(system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB")))
    expect_true(file.exists(system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam.bai", package = "PICB")))
})

test_that("PICBload loads BAM file with default parameters correctly", {
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    )
    expect_type(result, "list")
    expect_named(result, c("unique", "multi.primary", "multi.secondary"))
    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))
    expect_equal(length(result$unique), 31165)
    expect_equal(length(result$multi.primary), 47741)
    expect_equal(length(result$multi.secondary), 369219)
    expect_equal(ncol(mcols(result$unique)), 3)
    expect_equal(ncol(mcols(result$multi.primary)), 3)
    expect_equal(ncol(mcols(result$multi.secondary)), 3)
    expect_equal(min(start(result$unique)), 20063379)
    expect_equal(min(start(result$multi.primary)), 20006039)
    expect_equal(min(start(result$multi.secondary)), 20005922)
    expect_equal(max(end(result$unique)), 20984760)
    expect_equal(max(end(result$multi.primary)), 20980831)
    expect_equal(max(end(result$multi.secondary)), 20981429)
})

test_that("PICBload throws error when BAMFILE is NULL", {
    expect_error(
        PICBload(
        BAMFILE = NULL,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6"
        ),
        "Please provide full path to a .bam file!"
    )
})

test_that("PICBload throws error when REFERENCE.GENOME is NULL", {
    expect_error(
        PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = NULL
        ),
        "Please provide REFERENCE.GENOME"
    )
})

test_that("PICBload throws error when REFERENCE.GENOME is NULL", {
    expect_error(
        PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = NULL
        ),
        "Please provide REFERENCE.GENOME"
    )
})

test_that("PICBload handles SIMPLE.CIGAR = FALSE correctly", {
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        SIMPLE.CIGAR = FALSE
    )
    expect_type(result, "list")
    expect_named(result, c("unique", "multi.primary", "multi.secondary"))
    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))
    expect_equal(length(result$unique), 31207)
    expect_equal(length(result$multi.primary), 47844)
    expect_equal(length(result$multi.secondary), 369964)
    expect_equal(ncol(mcols(result$unique)), 3)
    expect_equal(ncol(mcols(result$multi.primary)), 3)
    expect_equal(ncol(mcols(result$multi.secondary)), 3)
})

test_that("PICBload loads secondary alignments when IS.SECONDARY.ALIGNMENT = TRUE", {
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        IS.SECONDARY.ALIGNMENT = TRUE
    )
    expect_type(result, "list")
    expect_named(result, c("multi.secondary"))
    expect_null(result$unique)
    expect_null(result$multi.primary)
    expect_true(inherits(result$multi.secondary, "GRanges"))
    expect_equal(length(result$multi.secondary), 369219)
    expect_equal(ncol(mcols(result$multi.secondary)), 3)
})

test_that("PICBload loads primary alignments when IS.SECONDARY.ALIGNMENT = FALSE", {
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        IS.SECONDARY.ALIGNMENT = FALSE
    )
    expect_type(result, "list")
    expect_named(result, c("unique", "multi.primary"))
    expect_null(result$multi.secondary)
    expect_true(inherits(result$unique, "GRanges"))
    expect_true(inherits(result$multi.primary, "GRanges"))
    expect_equal(length(result$unique), 31165)
    expect_equal(length(result$multi.primary), 47741)
    expect_equal(ncol(mcols(result$unique)), 3)
    expect_equal(ncol(mcols(result$multi.primary)), 3)
})

test_that("PICBload loads BAM file with default parameters the same with using the same reference genome but different methods", {
    expectedResult <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6"
    )
    # Using chromosome names and lengths
    seqnames <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
    seqlengths <- c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352)
    myGenome <- Seqinfo(seqnames = seqnames, seqlengths = seqlengths)

    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = myGenome
    )
    expect_equal(result, expectedResult)

    # Using genome name
    myGenome <- GenomeInfoDb::Seqinfo(genome = "dm6")
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = myGenome
    )
    expect_equal(result, expectedResult)

    # Using a FASTA file
    chr2L_seq <- BSgenome.Dmelanogaster.UCSC.dm6[["chr2L"]]
    chr2L_seq_set <- DNAStringSet(chr2L_seq)
    names(chr2L_seq_set) <- "chr2L"
    temp_fasta <- tempfile(fileext = ".fasta")
    writeXStringSet(chr2L_seq_set, temp_fasta)
    myGenome <- PICBloadfasta(FASTA.NAME = temp_fasta)
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = myGenome
    )
    expect_equal(result, expectedResult)
    unlink(temp_fasta)
})

test_that("PICBload filters alignments based on custom READ.SIZE.RANGE", {
    custom_size_range <- c(24, 27) # not recommended size range! Just here for testing
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        USE.SIZE.FILTER = TRUE,
        READ.SIZE.RANGE = custom_size_range
    )
    expect_type(result, "list")
    expect_named(result, c("unique", "multi.primary", "multi.secondary"))
    expect_true(all(sapply(result, function(x) inherits(x, "GRanges"))))
    expect_equal(length(result$unique), 28742)
    expect_equal(length(result$multi.primary), 44773)
    expect_equal(length(result$multi.secondary), 347254)
    expect_equal(ncol(mcols(result$unique)), 3)
    expect_equal(ncol(mcols(result$multi.primary)), 3)
    expect_equal(ncol(mcols(result$multi.secondary)), 3)
    # Check that all alignments have width within the specified range
    all_widths <- c(width(result$unique), width(result$multi.primary), width(result$multi.secondary))
    expect_true(all(all_widths >= custom_size_range[1] & all_widths <= custom_size_range[2]))
})

test_that("PICBload retrieves original sequences when GET.ORIGINAL.SEQUENCE = TRUE", {
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        GET.ORIGINAL.SEQUENCE = TRUE
    )
    expect_type(result, "list")
    expect_equal(length(result$unique), 31165)
    expect_equal(length(result$multi.primary), 47741)
    expect_equal(length(result$multi.secondary), 369219)
    expect_equal(ncol(mcols(result$unique)), 4)
    expect_equal(ncol(mcols(result$multi.primary)), 4)
    expect_equal(ncol(mcols(result$multi.secondary)), 4)
    # Check that 'seq' column is present in metadata
    expect_true("seq" %in% colnames(mcols(result$unique)))
    expect_true("seq" %in% colnames(mcols(result$multi.primary)))
    expect_true("seq" %in% colnames(mcols(result$multi.secondary)))
})

test_that("PICBload changes seqlevels style when SEQ.LEVELS.STYLE is specified", {
    custom_style <- "NCBI"
    result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        SEQ.LEVELS.STYLE = custom_style
    )
    expect_type(result, "list")
    # Check that seqlevels style is updated
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$unique), custom_style)
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$multi.primary), custom_style)
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$multi.secondary), custom_style)
})

test_that("PICBload works correctly with a combination of custom parameters", {
     result <- PICBload(
        BAMFILE = test_bam,
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        SIMPLE.CIGAR = FALSE,
        IS.SECONDARY.ALIGNMENT = FALSE,
        STANDARD.CONTIGS.ONLY = TRUE,
        PERFECT.MATCH.ONLY = TRUE,
        FILTER.BY.FLAG = FALSE,
        SELECT.FLAG = c(0, 16),
        USE.SIZE.FILTER = TRUE,
        READ.SIZE.RANGE = c(20, 30),
        TAGS = c("NH", "NM"),
        WHAT = c("flag", "qwidth"),
        SEQ.LEVELS.STYLE = "NCBI",
        GET.ORIGINAL.SEQUENCE = TRUE
    )
    expect_type(result, "list")
    expect_named(result, c("unique", "multi.primary"))
    # Check GRanges content
    expect_true(all(mcols(result$unique)$NH == 1))
    expect_true(all(mcols(result$multi.primary)$NH > 1))
    # Check 'seq' column
    expect_true("seq" %in% colnames(mcols(result$unique)))
    expect_true("seq" %in% colnames(mcols(result$multi.primary)))
    # Check read sizes
    all_widths_unique <- width(result$unique)
    all_widths_multi <- width(result$multi.primary)
    expect_true(all(all_widths_unique >= 20 & all_widths_unique <= 30))
    expect_true(all(all_widths_multi >= 20 & all_widths_multi <= 30))
    # Check that metadata columns are extended (seq and qwidth)
    expect_equal(ncol(mcols(result$unique)), 5)
    expect_equal(ncol(mcols(result$multi.primary)), 5)
    # Check seqlevels style
    expect_equal(GenomeInfoDb::seqlevelsStyle(result$unique), "NCBI")
})

test_that("PICBload VERBOSE parameter works correctly", {
    # Capture messages when VERBOSE is TRUE
    verbose_messages <- capture_messages({
    PICBload(
        BAMFILE = system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB"),
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSE = TRUE
    )
    })

    # Capture messages when VERBOSE is FALSE
    silent_messages <- capture_messages({
    PICBload(
        BAMFILE = system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB"),
        REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
        VERBOSE = FALSE
    )
    })

    # Test expectations
    expect_true(length(verbose_messages) > 0, "VERBOSE=TRUE should produce messages")
    expect_true(any(grepl("Processing", verbose_messages)), "Should contain 'Processing' message when VERBOSE=TRUE")
    expect_identical(length(silent_messages), 0L, "VERBOSE=FALSE should not produce any messages")
})
