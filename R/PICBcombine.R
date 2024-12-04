#' Combine Granges from few replicate into a set of repeatable intervals
#'
#' @param rangesList list of few PICBbuild outputs (with seeds, cores and clusters)
#' @param alignmentsList named list of PICBload outputs. Names correspond to replicate names.
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param TYPE.OF.REGION "clusters" by default
#' @param SEQ.LEVELS.STYLE style of chromosome names for BSgenome. "UCSC" by default.
#'
#' @return Granges object annotated for every replicate.
#' @export
#' @author Aleksandr Friman
PICBcombine <- function(rangesList = NULL, alignmentsList = NULL,
                        REFERENCE.GENOME = NULL, TYPE.OF.REGION = allalignments, SEQ.LEVELS.STYLE = "UCSC") {
    if (is.null(rangesList)) stop("Please provide rangesList !")
    if (is.null(alignmentsList)) stop("Please provide alignmentsList !")
    if (is.null(REFERENCE.GENOME)) stop("Please provide REFERENCE.GENOME !")
    # finding overlaps
    commonRanges <- rangesList[[1]][[TYPE.OF.REGION]]
    if (length(rangesList) > 1) {
        for (i in 2:length(rangesList)) {
            commonRanges <- GenomicRanges::intersect(commonRanges, rangesList[[i]][[TYPE.OF.REGION]])
        }
    }
    GenomicRanges::mcols(commonRanges) <- NULL
    # annotating mean coverage

    if (length(alignmentsList) > 0) {
        for (i in seq_along(alignmentsList)) {
            sampleName <- names(alignmentsList)[i]
            if (is.null(sampleName)) {
                sampleName <- i
            } else if (sampleName == "") {
                sampleName <- i
            }

            # names duplication check
            if (paste0("meanCoverage", sampleName) %in% colnames(GenomicRanges::mcols(commonRanges))) {
                stop("Sample names duplication! Check names of alignmentsList.")
            }
            if (TYPE.OF.REGION == allalignments) { # checking if type of cluster is relevant here
                    tmpRanges <- rangesList[[i]]
                    tmpRanges[[TYPE.OF.REGION]] <- commonRanges
                    tmpRanges <- PICBannotate(tmpRanges, alignmentsList[[i]],
                    REFERENCE.GENOME = REFERENCE.GENOME, REPLICATE.NAME = sampleName
                )
                commonRanges <- tmpRanges[[TYPE.OF.REGION]]
            } else {
                commonRanges <- PICBannotate(commonRanges, alignmentsList[[i]],
                    REFERENCE.GENOME = REFERENCE.GENOME, REPLICATE.NAME = sampleName
                )
            }
        }
    }

    commonRanges$hit <- NULL
    return(commonRanges)
}
