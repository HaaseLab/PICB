#' Export clustering object into an Excel file
#'
#' @param IN.RANGES clustering object to export
#' @param EXCEL.FILE.NAME file name to save
#'
#' @author Aleksandr Friman
#' @return no values returned
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
#' PICBexporttoexcel(
#'   IN.RANGES=outputOfPICBbuild, 
#'   EXCEL.FILE.NAME = "inst/extdata/myClusters_demonstration.xlsx")
PICBexporttoexcel <- function(IN.RANGES = NULL, EXCEL.FILE.NAME = NULL){
  allowedSheets = names(IN.RANGES)[ names(IN.RANGES) %in% c(uniqueonly, uniqueandprimary, allalignments) ]
  openxlsx::write.xlsx(IN.RANGES[ allowedSheets ], file=EXCEL.FILE.NAME)
}
