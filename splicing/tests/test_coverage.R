
context("Coverage to number of reads and back")

test_that("coverage to number of reads and back works", {

  library(splicing)

  gg <- createGene(list(c(1,500), c(1001,1100), c(951,1100), c(1301, 1400)),
                   list(c(1,2,4), c(3,4)))

  expr <- c(1/3, 2/3)
  rl <- 50
  reads <- simulateReads(gg, expression=expr, noReads=1000, readLength=rl)

  cv <- getCoverage(gg, reads=reads, exp=expr)
  nr <- coverage_to_no_reads(gg, exp=expr, coverage=cv, readLength=rl)
  expect_that(nr, equals(noReads(reads)))
  
})
