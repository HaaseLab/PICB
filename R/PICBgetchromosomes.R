#' Get SeqInfo object from standard non-circular chromosome names from your genome
#'
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6", or directly a SeqInfo object
#' @param SEQ.LEVELS.STYLE naming of chromosomes style. "UCSC" by default
#'
#' @return SeqInfo object with standard non-circular chromosome names
#' @export
#'
#' @author Aleksandr Friman
#' @author Franziska Ahrend
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' mySI <- PICBgetchromosomes("BSgenome.Dmelanogaster.UCSC.dm6", "UCSC")
PICBgetchromosomes <- function(REFERENCE.GENOME, SEQ.LEVELS.STYLE = DEFAULT.SEQ.LEVELS.STYLE) {
    if (typeof(REFERENCE.GENOME) == "character") {
        SI <- GenomeInfoDb::keepStandardChromosomes(GenomeInfoDb::seqinfo(x = BSgenome::getBSgenome(REFERENCE.GENOME)))
    } else {
        SI <- GenomeInfoDb::keepStandardChromosomes(REFERENCE.GENOME)
    }
    GenomeInfoDb::seqlevelsStyle(SI) <- SEQ.LEVELS.STYLE
    SI <- GenomeInfoDb::dropSeqlevels(x = SI, value = names(which(GenomeInfoDb::isCircular(SI))))
    return(SI)
}
