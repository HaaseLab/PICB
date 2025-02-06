#' Load and preprocess BAM files containing piRNA alignments
#'
#' @param BAMFILE name of the bam file to load. Should be sorted and indexed.
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param SIMPLE.CIGAR simpleCigar parameter of Rsamtools::ScanBamParam
#' @param IS.SECONDARY.ALIGNMENT defines loading of primary/secondary alignments. Default value NA loads both primary and secondary.
#' @param STANDARD.CONTIGS.ONLY use only standard chromosomes
#' @param PERFECT.MATCH.ONLY load only alignments without mismatches
#' @param FILTER.BY.FLAG enables filtering by flag. TRUE by default.
#' @param SELECT.FLAG vector of flags to use. Default value c(0,16, 272, 256).
#' @param USE.SIZE.FILTER enables filter by alignment size. True by default.
#' @param READ.SIZE.RANGE allowed alignment sizes. c(18,50) by default.
#' @param TAGS tags to import from bam file. c("NH","NM") by default.
#' @param WHAT "what" parameter of Rsamtools::ScanBamParam. c("flag") by default.
#' @param SEQ.LEVELS.STYLE naming of chromosomes style. "UCSC" by default.
#' @param GET.ORIGINAL.SEQUENCE adds "seq" to WHAT. False by default.
#' @param VERBOSE enables progress output. True by default.
#'
#' @author Pavol Genzor
#' @author Daniel Stoyko
#' @author Aleksandr Friman
#' @author Franziska Ahrend
#' @return list of GRanges objects named "unique" for unique mapping alignments,
#' "multi.primary" for primary multimapping alignments,
#' "multi.secondary" for secondary multimapping alignments
#' @export
#'
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' PICBload(
#'     BAMFILE = system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB"),
#'     REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#'     VERBOSE = FALSE
#' )
#'
PICBload <- function(
    ## INPUTS
    BAMFILE = NULL,
    REFERENCE.GENOME = NULL,
    ## OPTIONS
    SIMPLE.CIGAR = TRUE,
    IS.SECONDARY.ALIGNMENT = NA,
    STANDARD.CONTIGS.ONLY = TRUE,
    PERFECT.MATCH.ONLY = FALSE,
    FILTER.BY.FLAG = TRUE,
    SELECT.FLAG = c(0, 16, 272, 256),
    USE.SIZE.FILTER = TRUE,
    READ.SIZE.RANGE = c(18, 50),
    TAGS = c("NH", "NM"),
    WHAT = c("flag"),
    ## EXTRA OPTIONS
    SEQ.LEVELS.STYLE = "UCSC",
    GET.ORIGINAL.SEQUENCE = FALSE,
    VERBOSE = TRUE) {
    ## NOTE ON BAM INFO
    ## Flag: 256 = not primary alignment; 272 = reverse strand not primary alignment;
    ## Flag: 0 = forward unpaired unique alignment; 16 = reverse unpaired unique alignment
    ## Tags: NH:i:1 = unique alignment; NM = edit distance to the reference

    outputAlignments <- list()
    chrom_mismatch_flag <- FALSE
    ## check input
    if (is.null(BAMFILE)) stop("Please provide full path to a .bam file!")
    if (is.null(REFERENCE.GENOME)) stop("Please provide REFERENCE.GENOME")
    if (isTRUE(GET.ORIGINAL.SEQUENCE)) {
        WHAT <- c(WHAT, "seq")
    }
    ## for report

    ##
    if (VERBOSE) message("PICB v", utils::packageVersion("PICB"), " Processing ... ")
    ##

    justPrimaryOrSecondary <- function(IS.SECONDARY.ALIGNMENT) {
        ## PARAMETERS FOR LOADING BAM FILE
        if (isTRUE(STANDARD.CONTIGS.ONLY) && !is.na(SEQ.LEVELS.STYLE)) {
            genome_result <- PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)
            SI <- genome_result$SeqInfo
            chrom_mismatch_flag <- genome_result$chrom_mismatch

            REG.CHR <- GenomeInfoDb::seqnames(SI)
            BAM.FILE.HEADER <- Rsamtools::BamFile(BAMFILE)
            BAM.FILE.CHR <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BAM.FILE.HEADER))
            REG.CHR <- REG.CHR[REG.CHR %in% BAM.FILE.CHR]
            WHICH <- GenomicRanges::GRanges(
                seqnames = REG.CHR,
                ranges = IRanges::IRanges(start = rep(1, length(REG.CHR)), end = GenomeInfoDb::seqlengths(SI)[REG.CHR]),
                strand = rep("*", length(REG.CHR))
            )
            PARAM <- Rsamtools::ScanBamParam(
                flag = Rsamtools::scanBamFlag(
                    isUnmappedQuery = FALSE,
                    isSecondaryAlignment = IS.SECONDARY.ALIGNMENT
                ),
                tag = TAGS, simpleCigar = SIMPLE.CIGAR, what = WHAT, which = WHICH
            )
        } else {
            SI <- PICBgetchromosomes(REFERENCE.GENOME, NA)$SeqInfo
            PARAM <- Rsamtools::ScanBamParam(
                flag = Rsamtools::scanBamFlag(
                    isUnmappedQuery = FALSE,
                    isSecondaryAlignment = IS.SECONDARY.ALIGNMENT
                ),
                tag = TAGS, simpleCigar = SIMPLE.CIGAR, what = WHAT
            )
        }


        if (VERBOSE) {
            message("\nprepared loading parameters")
            message("\tTAGS:\t", toString(TAGS))
            message("\tCIGAR:\t", ifelse(isTRUE(SIMPLE.CIGAR), "simple cigar", "all cigar"))
            message("\tWHAT:\t", toString(WHAT))
            message("loading .bam file into GAlignments")

            if (is.na(IS.SECONDARY.ALIGNMENT)) {
                message("\n******")
                message("SLOW - Loading all reads")
                message(" => to load unique and primary alignments set")
                message(" => IS.SECONDARY.ALIGNMENT=FALSE")
                message("******\n")
            } else if (!IS.SECONDARY.ALIGNMENT) {
                message("\n******")
                message("Loading unique and primary alignments !!!")
                message("******\n")
            } else if (IS.SECONDARY.ALIGNMENT) {
                message("\n******")
                message("Loading secondary alignments only !!!")
                message("******\n")
            }
        }

        GA <- GenomicAlignments::readGAlignments(
            file = BAMFILE,
            use.names = TRUE,
            param = PARAM
        )
        # checking the tags consistency
        for (tagcheck in TAGS) {
            if (any(is.na(GenomicRanges::mcols(GA)[[tagcheck]]))) {
                stop("Tag ", tagcheck, " contains NA values. Check your bam file.")
            }
        }
        ## ***
        GA.IN <- length(GA)
        if (VERBOSE) message("\tIMPORTED: ", GA.IN)

        if (isTRUE(USE.SIZE.FILTER)) {
            if (VERBOSE) {
                message("\nfiltering by read size")
                message(
                "\tRANGE:\t",
                min(READ.SIZE.RANGE), "-", max(READ.SIZE.RANGE)
                )
            }
            GA <- GA[GenomicAlignments::width(GA) %in% seq(min(READ.SIZE.RANGE), max(READ.SIZE.RANGE), by = 1)]

            ## ***
            REMAINDER <- (length(GA) / GA.IN) * 100
            if (VERBOSE) message("\tREMAINDER: ", round(REMAINDER, digits = 2), " %")
        }



        if (isTRUE(FILTER.BY.FLAG)) {
            if (VERBOSE) message("\nfiltering based on flags")
            GARP <- GA[GenomicRanges::mcols(GA)[["flag"]] %in% SELECT.FLAG]
            # mcols(GARP)[["flag"]] <- NULL

            ## ***
            REMAINDER <- (length(GARP) / GA.IN) * 100
            if (VERBOSE) message("\tREMAINDER: ", round(REMAINDER, digits = 2), " %")
        } else {
            GARP <- GA
        }

        if (isTRUE(PERFECT.MATCH.ONLY)) {
            if (VERBOSE) message("\nremoving reads with mismatches")
            GARP <- GARP[GenomicRanges::mcols(GARP)[["NM"]] %in% c(0)]

            ## ***
            REMAINDER <- (length(GARP) / GA.IN) * 100
            if (VERBOSE) message("\tREMAINDER: ", round(REMAINDER, digits = 2), " %")
        }

        if (isTRUE(GET.ORIGINAL.SEQUENCE)) {
            if (VERBOSE) message("\nretrieving original read sequences")
            BAMSEQ <- GenomicRanges::mcols(GARP)[["seq"]]
            ISONMINUS <- as.logical(GenomicAlignments::strand(GARP) == "-")
            BAMSEQ[ISONMINUS] <- Biostrings::reverseComplement(BAMSEQ[ISONMINUS])
            GenomicRanges::mcols(GARP)[["seq"]] <- BAMSEQ
        }

        if (VERBOSE) message("\nconverting to GRanges")
        GARP.GR <- GenomicRanges::granges(GARP, use.names = TRUE, use.mcols = TRUE)
        # if no chromosomes mismatch and SEQ.LEVELS.STYLE is not NA, set the seqlevels style
        if (!chrom_mismatch_flag && !is.na(SEQ.LEVELS.STYLE)) {
            if (VERBOSE) message("\nSeqlevels style set to: ", SEQ.LEVELS.STYLE)
            GenomeInfoDb::seqlevelsStyle(GARP.GR) <- SEQ.LEVELS.STYLE
        }

        if (!all(unique(GenomeInfoDb::seqnames(GARP.GR)) %in% GenomeInfoDb::seqnames(SI))) {
            chrom_mismatch_flag <- TRUE
            warning("Some or all chromosome names in BAM file do not match those in REFERENCE.GENOME.",
                "\nNot filtered by chromosome names. Ensure REFERENCE.GENOME is compatible with the BAMFILE.")
        }

        return(GARP.GR)
    }
    if (VERBOSE) message("\nSorting into uniquemappers vs multimappers and primary vs secondary alignments")

    if (is.na(IS.SECONDARY.ALIGNMENT)) { # justPrimaryOrSecondary
        if (VERBOSE) message("Loading primary only")
        PrimaryAlignments <- justPrimaryOrSecondary(IS.SECONDARY.ALIGNMENT = FALSE)
        outputAlignments[["unique"]] <- PrimaryAlignments[GenomicRanges::mcols(PrimaryAlignments)[["NH"]] == 1]
        outputAlignments[["multi.primary"]] <- PrimaryAlignments[(GenomicRanges::mcols(PrimaryAlignments)[["NH"]] > 1)]
        if (VERBOSE) message("Loading secondary only")
        outputAlignments[["multi.secondary"]] <- justPrimaryOrSecondary(IS.SECONDARY.ALIGNMENT = TRUE)
    } else if (IS.SECONDARY.ALIGNMENT) {
        outputAlignments[["unique"]] <- NULL
        outputAlignments[["multi.primary"]] <- NULL
        if (VERBOSE) message("Loading secondary only")
        outputAlignments[["multi.secondary"]] <- justPrimaryOrSecondary(IS.SECONDARY.ALIGNMENT = TRUE)
    } else if (!IS.SECONDARY.ALIGNMENT) {
        if (VERBOSE) message("Loading primary only")
        PrimaryAlignments <- justPrimaryOrSecondary(IS.SECONDARY.ALIGNMENT = FALSE)
        outputAlignments[["unique"]] <- PrimaryAlignments[GenomicRanges::mcols(PrimaryAlignments)[["NH"]] == 1]
        outputAlignments[["multi.primary"]] <- PrimaryAlignments[(GenomicRanges::mcols(PrimaryAlignments)[["NH"]] > 1)]
        outputAlignments[["multi.secondary"]] <- NULL
    }
    ## RETURN
    if (chrom_mismatch_flag && !is.na(SEQ.LEVELS.STYLE)) {
        warning("No chromosome names in REFERENCE.GENOME match GenomeInfoDb::genomeStyle().",
        "\nKeeping all chromosomes given in REFERENCE.GENOME. No change of SEQ.LEVELS.STYLE.",
        "\nIf desired, filter REFERENCE.GENOME's SeqInfo for chromosomes (e.g. keeping only standard chromosomes).", 
        "\nTo not attempt to filter chromosomes or change SEQ.LEVELS.STYLE, set SEQ.LEVELS.STYLE to NA.")
    }
    if (VERBOSE) {
        message("\nDone!")
        message("")
    }
    return(outputAlignments)
}
