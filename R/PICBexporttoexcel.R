#' Export clustering object into an Excel file
#'
#' @param RANGES clustering object to export
#' @param EXCEL.FILE.NAME file name to save
#'
#' @author Aleksandr Friman
#' @return no values returned
#' @export
#'
#' @examples PICBexporttoexcel(myclusteringObject, "~/piRNAclustersFromUFOsamples.Area51.xlsx")
PICBexporttoexcel <- function(RANGES = NULL, EXCEL.FILE.NAME = NULL){
  allowedSheets = names(RANGES)[ names(RANGES) %in% c(uniqueonly, uniqueandprimary, allalignments) ]
  openxlsx::write.xlsx(RANGES[ allowedSheets ], file=EXCEL.FILE.NAME)
}
