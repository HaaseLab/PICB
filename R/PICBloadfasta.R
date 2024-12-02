#' Get SeqInfo object from a fasta file
#'
#' @param FASTA.NAME path to the fasta file
#'
#' @return SeqInfo object with all chromosome names and lengths from the fasta file
#' @export
#' @author Aleksandr Friman
#'
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#'
#' # create temporary fasta file
#' chr2L_seq <- BSgenome.Dmelanogaster.UCSC.dm6[["chr2L"]]
#' chr2L_seq_set <- DNAStringSet(chr2L_seq)
#' names(chr2L_seq_set) <- "chr2L"
#' temp_fasta <- tempfile(fileext = ".fasta")
#' writeXStringSet(chr2L_seq_set, temp_fasta)
#'
#' myGenome <- PICBloadfasta(FASTA.NAME = temp_fasta)
#' unlink(temp_fasta)
PICBloadfasta <- function(FASTA.NAME = NULL) {
  if (is.null(FASTA.NAME)) stop("Please provide FASTA.NAME !")
  FAdata <- seqinr::read.fasta(FASTA.NAME)
  FAnames <- names(FAdata)
  FAlengths <- lengths(FAdata)
  return(GenomeInfoDb::Seqinfo(seqnames = FAnames, seqlengths = FAlengths))
}
