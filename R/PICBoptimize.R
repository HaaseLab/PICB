#' Runs PICBbuild multiple times with provided parameters and returns optimization data frame.
#'
#' @param IN.ALIGNMENTS list of alignments from PICBload
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param LIBRARY.SIZE number of reads in the library. By default computed as number of unique mapping alignments + number of primary multimapping alignments.
#' @param VERBOSITY verbosity level 0/1/2/3. 2 by default.
#' @param PROVIDE.INFO.SEEDS.AND.CORES FALSE by default.
#' @param SEQ.LEVELS.STYLE naming of chromosomes style. "UCSC" by default.
#' @param ... rest of the parameters used by PICBbuild and provided as iteratable vectors
#'
#' @author Aleksandr Friman
#' @return Optimization values dataframe
#' @export
#'
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' myAlignmentsFromPICBload <- PICBload(
#'     BAMFILE = system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB"),
#'     REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#'     VERBOSE = FALSE
#' )
#'
#' PICBoptimize(
#'     IN.ALIGNMENTS = myAlignmentsFromPICBload,
#'     REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#'     MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2, 3, 4, 5)
#' )
PICBoptimize <- function(
    IN.ALIGNMENTS,
    REFERENCE.GENOME,
    LIBRARY.SIZE = length(IN.ALIGNMENTS$unique) + length(IN.ALIGNMENTS$multi.primary),
    VERBOSITY = 2,
    PROVIDE.INFO.SEEDS.AND.CORES = FALSE,
    SEQ.LEVELS.STYLE = "UCSC",
    ...) {

    if (VERBOSITY > 0) message("PICB v", utils::packageVersion("PICB"), " Starting ... ")

    # checking the inputs
    if (VERBOSITY > 0) {
        message("Checking the inputs")
    }
    if (is.null(IN.ALIGNMENTS)) stop("Please provide IN.ALIGNMENTS !")
    if (is.null(REFERENCE.GENOME)) stop("Please provide REFERENCE.GENOME !")

    numberOfallReadsExplained <- function(gr, alignments) {
        foundnames <- c()
        if (typeof(alignments) == "list") {
            for (clmn in c("unique", "multi.primary", "multi.secondary")) {
                usedAlignments <- IRanges::subsetByOverlaps(alignments[[clmn]], gr)
                foundnames <- c(foundnames, names(usedAlignments))
            }
        } else {
            usedAlignments <- IRanges::subsetByOverlaps(alignments, gr)
            foundnames <- c(foundnames, names(usedAlignments))
        }
        numbr <- length(unique(foundnames))
        return(numbr)
    }

    inputArgs <- as.list(match.call(expand.dots = TRUE))
    PICBbuildArgs <- formals(PICBbuild)
    ArgsToIterate <- NULL
    for (argName in names(PICBbuildArgs)) {
        # checking if the argument need to be iterated over
        if (is.numeric(PICBbuildArgs[[argName]]) ||
            (is.language(PICBbuildArgs[[argName]]) && PICBbuildArgs[[argName]] != "")) {
            # iteration is possible
            if (argName %in% names(inputArgs)) {
                # need to iterate
                if (length(inputArgs[[argName]]) > 1) {
                    if (VERBOSITY > 0) {
                        message("Adding ", argName, " to iteration targets")
                    }
                    if (is.null(ArgsToIterate)) {
                        ArgsToIterate <- argName
                    } else {
                        ArgsToIterate <- c(ArgsToIterate, argName)
                    }
                }
            }
        }
    }

    
    if (is.null(ArgsToIterate)) {
        stop("Provide arguments to iterave over. See example.")
    }

    totalReads <- LIBRARY.SIZE
    if (length(unique(c(
        names(IN.ALIGNMENTS$unique),
        names(IN.ALIGNMENTS$multi.primary)
    ))) != totalReads) {
        warning("The total number of primary alignments is not equal to the total number of read names. This discrepancy may occur when using a subset of data, such as in the PICB demonstration, where you specified the LIBRARY.SIZE but the actual number of reads is reduced. For full piRNA datasets, this indicates a potential issue with read assignment and your numbers of explained reads may be wrong.")
    }
    # building the parameter combinations
    if (VERBOSITY > 0) {
        message("Building the parameter combinations")
    }
    listOfArgLists <- list()
    for (argName in ArgsToIterate) {
        firstArg <- length(listOfArgLists) == 0
        argVals <- eval(inputArgs[[argName]])
        if (firstArg) { # first arg
            for (curVal in argVals) {
                tmpList <- list()
                tmpList[[argName]] <- curVal
                listOfArgLists[[length(listOfArgLists) + 1]] <- tmpList
            }
        } else {
            # add the argName to each combination
            newlistOfArgLists <- list()
            for (combList in listOfArgLists) {
                for (curVal in argVals) {
                    combList[[argName]] <- curVal
                    newlistOfArgLists[[length(newlistOfArgLists) + 1]] <- combList
                }
            }
            listOfArgLists <- newlistOfArgLists
        }
    }
    IN.ALIGNMENTS.ALL <- c(
        IN.ALIGNMENTS$unique, IN.ALIGNMENTS$multi.primary,
        IN.ALIGNMENTS$multi.secondary
    )
    # preparing output dataframe
    outDF <- data.frame(combinationNumber = seq_along(listOfArgLists))
    for (argName in names(listOfArgLists[[1]])) {
        # adding columns we iterate over
        outDF[[argName]] <- NA
    }
    if (PROVIDE.INFO.SEEDS.AND.CORES) {
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("number.of.", locustype)]] <- NA
        }
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("total.width.", locustype)]] <- NA
        }
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("reads.explained.by.", locustype)]] <- NA
        }
    } else {
        outDF[[paste0("number.of.clusters")]] <- NA
        outDF[[paste0("total.width.clusters")]] <- NA
        outDF[[paste0("reads.explained.by.clusters")]] <- NA
    }
    # running PICBbuild for every possible combation
    i <- 1
    for (combList in listOfArgLists) {
        if (VERBOSITY > 0) {
            message("Iteration ", i, " out of ", length(listOfArgLists))
        }
        for (argName in names(combList)) {
            # adding values we iterate over
            outDF[[argName]][i] <- combList[[argName]]
        }
        listToRun <- combList
        listToRun[["IN.ALIGNMENTS"]] <- IN.ALIGNMENTS
        listToRun[["REFERENCE.GENOME"]] <- REFERENCE.GENOME
        listToRun[["VERBOSITY"]] <- VERBOSITY
        if (VERBOSITY > 1) {
            message("Running PICBbuild")
        }
        tmrLoci <- do.call(what = PICBbuild, args = listToRun)
        if (VERBOSITY > 1) {
            message("PICBbuild done")
        }
        if (PROVIDE.INFO.SEEDS.AND.CORES) {
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("number.of.", locustype)]][i] <- length(tmrLoci[[locustype]])
            outDF[[paste0("total.width.", locustype)]][i] <- as.numeric(sum(GenomicRanges::width(tmrLoci[[locustype]])))
            outDF[[paste0("reads.explained.by.", locustype)]][i] <- as.numeric(numberOfallReadsExplained(tmrLoci[[locustype]], IN.ALIGNMENTS.ALL))
        }
        } else {
            outDF[["number.of.clusters"]][i] <- length(tmrLoci[[allalignments]])
            outDF[["total.width.clusters"]][i] <- as.numeric(sum(GenomicRanges::width(tmrLoci[[allalignments]])))
            outDF[["reads.explained.by.clusters"]][i] <- as.numeric(numberOfallReadsExplained(tmrLoci[[allalignments]], IN.ALIGNMENTS.ALL))
        }
        i <- i + 1
    }
    outDF$combinationNumber <- NULL
    # adding normaized values
    if (VERBOSITY > 0) {
        message("Values normalization")
    }
    if (PROVIDE.INFO.SEEDS.AND.CORES) {
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("fraction.of.library.explained.by.", locustype)]] <- outDF[[paste0("reads.explained.by.", locustype)]] / totalReads
        }
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
            outDF[[paste0("mean.RPKM.", locustype)]] <- outDF[[paste0("reads.explained.by.", locustype)]] * 1e9 / (outDF[[paste0("total.width.", locustype)]] * totalReads)
            outDF[[paste0("mean.RPKM.", locustype)]][is.na(outDF[[paste0("mean.RPKM.", locustype)]])] <- 0
        }
        SI <- PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)$SeqInfo
        GenomeSize <- sum(GenomeInfoDb::seqlengths(SI))
        for (locustype in c(uniqueonly, uniqueandprimary, allalignments)) {
        outDF[[paste0("fraction.of.genome.space.", locustype)]] <- outDF[[paste0("total.width.", locustype)]] / (2 * GenomeSize)
        }
    } else {
        outDF[["fraction.of.library.explained.by.clusters"]] <- outDF[["reads.explained.by.clusters"]] / totalReads
        outDF[["mean.RPKM.clusters"]] <- outDF[["reads.explained.by.clusters"]] * 1e9 / (outDF[["total.width.clusters"]] * totalReads)
        outDF[["mean.RPKM.clusters"]][is.na(outDF[["mean.RPKM.clusters"]])] <- 0
        SI <- PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)$SeqInfo
        GenomeSize <- sum(GenomeInfoDb::seqlengths(SI))
        outDF[["fraction.of.genome.space.clusters"]] <- outDF[["total.width.clusters"]] / (2 * GenomeSize)
    }
    return(outDF)
}