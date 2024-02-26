#' Get SeqInfo object from a fasta file
#'
#' @param FASTA.NAME path to the fasta file
#'
#' @return SeqInfo object with all chromosome names and lengths from the fasta file
#' @export
#'
#' @examples mySI<-PICBloadfasta("~/path/to/your.fasta")
PICBloadfasta<-function(FASTA.NAME = NULL){
  if(is.null(FASTA.NAME)) stop("Please provide FASTA.NAME !")
  FAdata<-seqinr::read.fasta(FASTA.NAME)
  FAnames<-names(FAdata)
  FAlengths<-lengths(FAdata)
  return(GenomeInfoDb::Seqinfo(seqnames = FAnames, seqlengths = FAlengths))
}

