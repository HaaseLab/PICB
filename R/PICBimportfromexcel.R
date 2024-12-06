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
PICBimportfromexcel <- function(EXCEL.FILE.NAME = NULL, REFERENCE.GENOME = NULL,
        SEQ.LEVELS.STYLE = DEFAULT.SEQ.LEVELS.STYLE) {
    availbsheets <- openxlsx::getSheetNames(EXCEL.FILE.NAME)
    output <- list()
    sheetsToCheck <- intersect(availbsheets, c(uniqueonly, uniqueandprimary, allalignments))
    for (t in sheetsToCheck) {
        output[[t]] <- GenomicRanges::GRanges(openxlsx::read.xlsx(EXCEL.FILE.NAME, sheet = t))
        
    }
    if (is.null(REFERENCE.GENOME)){
        warning("REFERENCE.GENOME not providen. Cannot add chromosome length information.")
        return(output)
    } 
    SI <- PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)
    for (t in names(output)) {
        GenomicRanges::seqinfo(output[[t]]) <- SI
    }
    return(output)
}
