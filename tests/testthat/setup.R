# tests/testthat/setup.R
library("BSgenome.Dmelanogaster.UCSC.dm6")

# Common test data setup
test_bam <- system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB")

# Load alignments once
test_alignments <- PICBload(
    BAMFILE = test_bam,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSE = FALSE
)

# Build ranges (seeds, cores, clusters) once
test_ranges <- PICBbuild(
    IN.ALIGNMENTS = test_alignments,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSITY = 0
)

test_optimization <- PICBoptimize(
    IN.ALIGNMENTS = test_alignments,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSITY = 0,
    MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2)
)