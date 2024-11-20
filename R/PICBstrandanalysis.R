#' Add sense/antisense ratio of piRNAs per piRNAcluster
#'
#' @param IN.ALIGNMENTS list of alignments from PICBload
#' @param IN.RANGES single GRanges object (seeds, cores or clusters from PICBbuild)
#' @param VERBOSE enables progress output. True by default.
#'
#' @return Granges object with an additional annotation column
#' @export
#'
#' @examples outputOfPICBbuild$clusters<-PICBstrandanalysis(outputOfPICBload, outputOfPICBbuild$clusters)
PICBstrandanalysis<-function(IN.ALIGNMENTS,
                             IN.RANGES,
                             VERBOSE=TRUE){
    
    if(is.null(IN.ALIGNMENTS)) stop("Please provide IN.ALIGNMENTS !")
    if(is.null(IN.RANGES)) stop("Please provide IN.RANGES !")
    if(is.list(IN.RANGES)) stop("IN.RANGES must be a single GRanges object (seeds, cores or clusters)!")


    if (VERBOSE==TRUE) message(paste("PICB v", packageVersion("PICB"),"Starting ... "))

    if (VERBOSE==TRUE) message("Using unique alignments only")
    # use unique alignments only
    IN.ALIGNMENTS <- IN.ALIGNMENTS$unique[IN.ALIGNMENTS$unique$NH == 1]
    # get strandedness
    IN.ALIGNMENTS.PLUS <- IN.ALIGNMENTS[GenomicRanges::strand(IN.ALIGNMENTS) == "+"]
    IN.ALIGNMENTS.MINUS <- IN.ALIGNMENTS[GenomicRanges::strand(IN.ALIGNMENTS) == "-"]

    if (VERBOSE==TRUE) message("Finding overlaps of clusters with piRNAs from plus strand")
    DT_plus <- as.data.frame(GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.PLUS, ignore.strand = TRUE))
    DT_plus.2 <- aggregate(subjectHits ~ queryHits, data = DT_plus, FUN = length)
    # Add numbers of plus-strand piRNAs to clusters
    IN.RANGES$plus_piRNAs <- 0
    IN.RANGES$plus_piRNAs[DT_plus.2$queryHits] <- DT_plus.2$subjectHits

    if (VERBOSE==TRUE) message("Finding overlaps of clusters with piRNAs from minus strand")
    DT_minus <- as.data.frame(GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.MINUS, ignore.strand = TRUE))
    DT_minus.2 <- aggregate(subjectHits ~ queryHits, data = DT_minus, FUN = length)
    
    # Add numbers of plus-strand piRNAs to clusters
    IN.RANGES$minus_piRNAs <- 0
    IN.RANGES$minus_piRNAs[DT_minus.2$queryHits] <- DT_minus.2$subjectHits

    #to avoid deviding by 0, add 1 piRNA to each cluster (plus and minus strands)
    IN.RANGES$plus_piRNAs <- IN.RANGES$plus_piRNAs + 1
    IN.RANGES$minus_piRNAs <- IN.RANGES$minus_piRNAs + 1

    if (VERBOSE==TRUE) message("Adding ratio of sense/antisense piRNAs per cluster")
    IN.RANGES$s_as_ratio <- 1
    IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$s_as_ratio <- (IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$plus_piRNAs/IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$minus_piRNAs)
    IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$s_as_ratio <- (IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$minus_piRNAs/IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$plus_piRNAs)

    #remove plus_piRNAs and minus_piRNAs
    IN.RANGES$plus_piRNAs <- NULL
    IN.RANGES$minus_piRNAs <- NULL

    if (VERBOSE==TRUE) message("Done!")

    return(IN.RANGES)
}


