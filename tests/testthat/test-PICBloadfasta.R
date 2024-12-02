# tests/testthat/test-PICBloadfasta.R
test_that("PICBloadfasta loads a valid FASTA file correctly", {
  temp_fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">chr1", "ATCG", ">chr2", "GGT"), temp_fasta)

  # Expected output (mocked for illustration)
  expected_output <- GenomeInfoDb::Seqinfo(
    seqnames = c("chr1", "chr2"),
    seqlengths = c(4, 3),
    isCircular = c(NA, NA),
    genome = c(NA, NA)
  )

  expect_equal(PICBloadfasta(temp_fasta), expected_output)
  unlink(temp_fasta)
})

test_that("PICBloadfasta handles invalid input of FASTA file", {
  expect_error(PICBloadfasta(NULL), "Please provide FASTA.NAME !")

  temp_fasta <- tempfile(fileext = ".fasta")
  writeLines("", temp_fasta)
  expect_error(PICBloadfasta(temp_fasta), "no line starting with a > character found")
  writeLines("ATCG", temp_fasta)
  expect_error(PICBloadfasta(temp_fasta), "no line starting with a > character found")
  unlink(temp_fasta)
})
