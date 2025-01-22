test_bam <- system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB")
library("BSgenome.Dmelanogaster.UCSC.dm6")
myAlignments <- PICBload(
    BAMFILE = test_bam,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSE = FALSE
)
myRanges <- PICBbuild(
    IN.ALIGNMENTS = myAlignments,
    REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6",
    VERBOSITY = 0
)

test_that("PICBexporttoexcel creates correct Excel files", {
  # Test full export with multiple sheets
  temp_file <- tempfile(fileext = ".xlsx")
  PICBexporttoexcel(IN.RANGES = myRanges, EXCEL.FILE.NAME = temp_file)
  
  # Read the Excel file back
  xlsx_data <- openxlsx::loadWorkbook(temp_file)
  
  # Check number of sheets
  expect_equal(length(names(xlsx_data)), 3)
  
  # Check content of sheets
  for (sheet in names(xlsx_data)) { 
    sheet_data <- openxlsx::read.xlsx(temp_file, sheet = sheet)
    expect_gt(nrow(sheet_data), 1) # More than just header row
    expect_gt(ncol(sheet_data), 0) # Has columns
  }
  
  # Test single GRanges export
  single_temp_file <- tempfile(fileext = ".xlsx")
  PICBexporttoexcel(IN.RANGES = myRanges$clusters, 
                    EXCEL.FILE.NAME = single_temp_file)
  
  # Read single sheet Excel file back
  single_xlsx_data <- openxlsx::loadWorkbook(single_temp_file)
  
  # Check single sheet export
  expect_equal(length(names(single_xlsx_data)), 1)
  expect_equal(names(single_xlsx_data), "clusters")
  
  sheet_data <- openxlsx::read.xlsx(single_temp_file, sheet = 1)
  expect_gt(nrow(sheet_data), 1)
  expect_gt(ncol(sheet_data), 0)
  
  # Clean up
  unlink(temp_file)
  unlink(single_temp_file)
})
