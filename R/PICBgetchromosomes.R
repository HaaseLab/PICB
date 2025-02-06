#' Get SeqInfo object from standard non-circular chromosome names from your genome
#'
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6", or directly a SeqInfo object.
#' @param SEQ.LEVELS.STYLE naming of chromosomes style. "UCSC" by default.
#'
#' @return SeqInfo object with standard non-circular chromosome names
#' @export
#'
#' @author Aleksandr Friman
#' @author Franziska Ahrend
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' mySI <- PICBgetchromosomes("BSgenome.Dmelanogaster.UCSC.dm6", "UCSC")
PICBgetchromosomes <- function(REFERENCE.GENOME, SEQ.LEVELS.STYLE = "UCSC") {    
    # Get SeqInfo object first, regardless of SEQ.LEVELS.STYLE
    if (methods::is(REFERENCE.GENOME, "character")) {  #BSgenome object
        SI <- GenomeInfoDb::seqinfo(x = eval(parse(text = REFERENCE.GENOME)))
    } else if (methods::is(REFERENCE.GENOME, "Seqinfo")) {
        SI <- REFERENCE.GENOME
    } else {
        stop("REFERENCE.GENOME must be a character BSgenome object (e.g. 'BSgenome.Dmelanogaster.UCSC.dm6') or a Seqinfo object")
    }
    
    if (is.na(SEQ.LEVELS.STYLE)) {
        # If SEQ.LEVELS.STYLE is NA, return full SeqInfo without filtering
        return(list(SeqInfo = SI, chrom_mismatch = TRUE))
    }

    # Filter for standard chromosomes
    SI <- GenomeInfoDb::keepStandardChromosomes(SI) 
    
    # Check if any standard chromosomes were found
    if (length(SI) == 0) {
        return(list(SeqInfo = REFERENCE.GENOME, chrom_mismatch = TRUE))
    }
    
    # Update style and drop circular chromosomes
    GenomeInfoDb::seqlevelsStyle(SI) <- SEQ.LEVELS.STYLE
    SI <- GenomeInfoDb::dropSeqlevels(x = SI, value = names(which(GenomeInfoDb::isCircular(SI))))

    return(list(SeqInfo = SI, chrom_mismatch = FALSE))
}
