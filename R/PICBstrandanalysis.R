#' Add sense/antisense ratio of unique piRNAs per piRNAcluster
#'
#' @param IN.ALIGNMENTS list of alignments from PICBload
#' @param IN.RANGES single GRanges object (seeds, cores or clusters from PICBbuild)
#' @param VERBOSE enables progress output. True by default.
#'
#' @return Granges object with an additional annotation column
#' @author Parthena Konstantinidou
#' @author Zuzana Loubalova
#' @author Franziska Ahrend
#' @export
#'
#' @examples 
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' myAlignmentsFromPICBload <- PICBload(
#'   BAMFILE=system.file("extdata", "Fly_Ov1_chr2L_20To21mb.bam", package = "PICB"), 
#'   REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6", 
#'   VERBOSE=FALSE
#' )
#' 
#' outputOfPICBbuild <- PICBbuild(
#'   IN.ALIGNMENTS = myAlignmentsFromPICBload,
#'   REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6", 
#'   VERBOSITY=0
#' )
#' 
#' outputOfPICBbuild$clusters<-PICBstrandanalysis(
#'   IN.ALIGNMENTS = myAlignmentsFromPICBload, 
#'   IN.RANGES = outputOfPICBbuild$clusters
#' )
PICBstrandanalysis<-function(IN.ALIGNMENTS,
                             IN.RANGES,
                             VERBOSE=TRUE){
    
    if(is.null(IN.ALIGNMENTS) || length(IN.ALIGNMENTS$unique) == 0) stop("IN.ALIGNMENTS is either NULL or contains no rows. Please provide valid IN.ALIGNMENTS from PICBload!")
    if(is.null(IN.RANGES) || length(IN.RANGES)==0) stop("IN.RANGES is either NULL or contains no rows. Use the GRanges object seeds, cores, or clusters from PICBbuild to create valid IN.RANGES!")
    if(is.list(IN.RANGES)) stop("IN.RANGES must be a single GRanges object: seeds, cores, or clusters from PICBbuild!")
    

    if (VERBOSE==TRUE) message("PICB v", utils::packageVersion("PICB")," Starting with strand analysis ... ")

    if (VERBOSE==TRUE) message("\tUsing unique alignments only")
    # use unique alignments only
    IN.ALIGNMENTS <- IN.ALIGNMENTS$unique[IN.ALIGNMENTS$unique$NH == 1]
    # get strandedness
    IN.ALIGNMENTS.PLUS <- IN.ALIGNMENTS[GenomicRanges::strand(IN.ALIGNMENTS) == "+"]
    IN.ALIGNMENTS.MINUS <- IN.ALIGNMENTS[GenomicRanges::strand(IN.ALIGNMENTS) == "-"]
    
    if (VERBOSE==TRUE) message("\tFinding overlaps of clusters with piRNAs from plus strand")
    #DT_plus <- as.data.frame(GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.PLUS, ignore.strand = TRUE))
    overlaps_plus <- GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.PLUS, ignore.strand = TRUE)
    DT_plus <- data.frame(
        queryHits = S4Vectors::queryHits(overlaps_plus),
        subjectHits = S4Vectors::subjectHits(overlaps_plus)
    )
    if (nrow(DT_plus) == 0) {
        warning("No overlaps detected between IN.ALIGNMENTS and clusters on the plus strand. Please verify that IN.ALIGNMENTS and IN.RANGES are valid and correctly correspond to each other. Continuing ...")
        DT_plus.2 <- data.frame(queryHits = integer(0), subjectHits = integer(0))
    } else {
        DT_plus.2 <- stats::aggregate(subjectHits ~ queryHits, data = DT_plus, FUN = length)
    }
    # Add numbers of plus-strand piRNAs to clusters
    IN.RANGES$plus_piRNAs <- 0
    IN.RANGES$plus_piRNAs[DT_plus.2$queryHits] <- DT_plus.2$subjectHits

    if (VERBOSE==TRUE) message("\tFinding overlaps of clusters with piRNAs from minus strand")
    #DT_minus <- as.data.frame(GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.MINUS, ignore.strand = TRUE))
    overlaps_minus <- GenomicRanges::findOverlaps(IN.RANGES, IN.ALIGNMENTS.MINUS, ignore.strand = TRUE)
    DT_minus <- data.frame(
        queryHits = S4Vectors::queryHits(overlaps_minus),
        subjectHits = S4Vectors::subjectHits(overlaps_minus)
    )
    if (nrow(DT_minus) == 0) {
        warning("No overlaps detected between IN.ALIGNMENTS and clusters on the minus strand. Please verify that IN.ALIGNMENTS and IN.RANGES are valid and correctly correspond to each other. Continuing ...")
        DT_minus.2 <- data.frame(queryHits = integer(0), subjectHits = integer(0))
    } else {
        DT_minus.2 <- stats::aggregate(subjectHits ~ queryHits, data = DT_minus, FUN = length)
    }
    # Add numbers of plus-strand piRNAs to clusters
    IN.RANGES$minus_piRNAs <- 0
    IN.RANGES$minus_piRNAs[DT_minus.2$queryHits] <- DT_minus.2$subjectHits

    #to avoid deviding by 0, add 1 piRNA to each cluster (plus and minus strands)
    IN.RANGES$plus_piRNAs <- IN.RANGES$plus_piRNAs + 1
    IN.RANGES$minus_piRNAs <- IN.RANGES$minus_piRNAs + 1

    if (VERBOSE==TRUE) message("\tAdding ratio of sense/antisense piRNAs per cluster")
    IN.RANGES$s_as_ratio <- 1
    IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$s_as_ratio <- (IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$plus_piRNAs/IN.RANGES[GenomicRanges::strand(IN.RANGES) == "+"]$minus_piRNAs)
    IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$s_as_ratio <- (IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$minus_piRNAs/IN.RANGES[GenomicRanges::strand(IN.RANGES) == "-"]$plus_piRNAs)

    #remove plus_piRNAs and minus_piRNAs
    IN.RANGES$plus_piRNAs <- NULL
    IN.RANGES$minus_piRNAs <- NULL

    if (VERBOSE==TRUE) message("Done!")

    return(IN.RANGES)
}


