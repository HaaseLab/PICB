test_that("PICBimportfromexcel correctly imports Excel files", {
  # Setup: Create test Excel files using export
  temp_file <- tempfile(fileext = ".xlsx")
  single_temp_file <- tempfile(fileext = ".xlsx")
  
  # Export test data
  PICBexporttoexcel(IN.RANGES = test_ranges, EXCEL.FILE.NAME = temp_file)
  PICBexporttoexcel(IN.RANGES = list(clusters = test_ranges$clusters), 
                    EXCEL.FILE.NAME = single_temp_file)
  
  # Test importing multi-sheet Excel
  imported_ranges <- PICBimportfromexcel(EXCEL.FILE.NAME = temp_file)
  
  # Check structure
  expect_type(imported_ranges, "list")
  expect_equal(length(imported_ranges), 3)
  expect_true(all(c("seeds", "cores", "clusters") %in% names(imported_ranges)))
  
  # Check that each element is a GRanges object
  for (name in names(imported_ranges)) {
    expect_s4_class(imported_ranges[[name]], "GRanges")
    # Compare with original
    expect_equal(length(imported_ranges[[name]]), length(test_ranges[[name]]))
    expect_equal(as.character(seqnames(imported_ranges[[name]])), as.character(seqnames(test_ranges[[name]])))
    expect_equal(ranges(imported_ranges[[name]]), ranges(test_ranges[[name]]))
  }
  
  # Test importing single-sheet Excel
  imported_single <- PICBimportfromexcel(EXCEL.FILE.NAME = single_temp_file)
  
  # Check structure of single import
  expect_type(imported_single, "list")
  expect_equal(length(imported_single), 1)
  expect_true("clusters" %in% names(imported_single))
  
  # Check content of single import
  expect_s4_class(imported_single$clusters, "GRanges")
  expect_equal(length(imported_single$clusters), length(test_ranges$clusters))
  expect_equal(as.character(seqnames(imported_single$clusters)), as.character(seqnames(test_ranges$clusters)))
  expect_equal(ranges(imported_single$clusters), ranges(test_ranges$clusters))
  
  # Clean up
  unlink(temp_file)
  unlink(single_temp_file)
})
