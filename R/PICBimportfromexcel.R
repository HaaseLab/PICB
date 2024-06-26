#' Import clustering object from an Excel file
#'
#' @param EXCEL.FILE.NAME file name to import from
#'
#' @author Aleksandr Friman
#'
#' @return list of annotated Granges objects named "seeds" for seeds,
#' "cores" for cores,
#' "clusters" for clusters
#' @export
#'
#' @examples PICBimportfromexcel(="~/piRNAclustersFromUFOsamples.Area51.xlsx")
PICBimportfromexcel <- function(EXCEL.FILE.NAME = NULL){
  availbsheets=openxlsx::getSheetNames(EXCEL.FILE.NAME)
  output=list()
  sheetsToCheck=intersect(availbsheets, c(uniqueonly, uniqueandprimary, allalignments))
  for (t in sheetsToCheck){
    output[[t]]<-GenomicRanges::GRanges(openxlsx::read.xlsx(EXCEL.FILE.NAME, sheet = t))
  }
  return(output)
}

