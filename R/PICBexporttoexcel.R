#' Export cluster object into an Excel file
#'
#' @param IN.RANGES clustering object to export
#' @param EXCEL.FILE.NAME file name to save
#'
#' @author Aleksandr Friman
#' @author Franziska Ahrend
#' @return no values returned
#' @export
#'
#' @examples
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' myAlignmentsFromPICBload <- PICBload(
#'     BAMFILE = system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB"),
#'     REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#'     VERBOSE = FALSE
#' )
#'
#' outputOfPICBbuild <- PICBbuild(
#'     IN.ALIGNMENTS = myAlignmentsFromPICBload,
#'     REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
#'     LIBRARY.SIZE = 12799826, #usually not necessary
#'     VERBOSITY = 0
#' )
#' PICBexporttoexcel(
#'     IN.RANGES = outputOfPICBbuild,
#'     EXCEL.FILE.NAME = "inst/extdata/myClusters_demonstration.xlsx"
#' )
PICBexporttoexcel <- function(IN.RANGES = NULL, EXCEL.FILE.NAME = NULL) {
    # Check if IN.RANGES is a single GRanges object and convert to list
    if (inherits(IN.RANGES, "GRanges")) {
        message("IN.RANGES is provided as a single GRanges object. Converting to a list with the name 'clusters'. If you intend to save seeds or cores, please input the data as a list (e.g., list(seeds = IN.RANGES)).")
        IN.RANGES <- list(clusters = IN.RANGES)
    }
    valid_names <- c(uniqueonly, uniqueandprimary, allalignments)
    invalid_names <- setdiff(names(IN.RANGES), valid_names) 
    if (length(invalid_names) > 0) {
        stop("Invalid names detected in IN.RANGES: ", paste(invalid_names, collapse = ", "), ". Expected names are 'seeds', 'cores', and/or 'clusters'.")
    }

    allowedSheets <- names(IN.RANGES)[names(IN.RANGES) %in% c(uniqueonly, uniqueandprimary, allalignments)]
    openxlsx::write.xlsx(IN.RANGES[allowedSheets], file = EXCEL.FILE.NAME)
    message("Saved to Excel.")
}
