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
  output=list()
  output$seeds<-GenomicRanges::GRanges(openxlsx::read.xlsx(EXCEL.FILE.NAME, sheet = uniqueonly))
  output$cores<-GenomicRanges::GRanges(openxlsx::read.xlsx(EXCEL.FILE.NAME, sheet = uniqueandprimary))
  output$clusters<-GenomicRanges::GRanges(openxlsx::read.xlsx(EXCEL.FILE.NAME, sheet = uniqueandprimary))
  return(output)
}

