#' Get SeqInfo object from BSgenome name
#'
#' @param BS.SPECIES name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param SEQ.LEVELS.STYLE naming of chromosomes style. "UCSC" by default
#'
#' @return SeqInfo object with standard non-circular chromosome names
#' @export
#'
#' @author Aleksandr Friman
#' @examples 
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' mySI<-PICBgetchromosomes("BSgenome.Dmelanogaster.UCSC.dm6", "UCSC")
PICBgetchromosomes<-function(BS.SPECIES, SEQ.LEVELS.STYLE = "UCSC"){
  #library(BS.SPECIES, character.only = TRUE)
  SI <- GenomeInfoDb::keepStandardChromosomes(GenomeInfoDb::seqinfo(x = eval(parse(text = BS.SPECIES))))
  GenomeInfoDb::seqlevelsStyle(SI) <- SEQ.LEVELS.STYLE
  SI <- GenomeInfoDb::dropSeqlevels(x = SI, value = names(which(GenomeInfoDb::isCircular(SI))))
  return(SI)
}
