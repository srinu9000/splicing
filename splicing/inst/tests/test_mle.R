
context("mle solver")

test_that("mle solver works", {
  library(splicing)
  set.seed(42)

  ## Skipped exon
  g <- createGene(list(c(1,100), c(201,300), c(401,500)),
                  list(c(1,2,3), c(1,3)))

  reads <- simulateReads(g, expression=c(1/3, 2/3), noReads=2000,
                         readLength=30)

  sol <- mleIso(g, reads=reads)

  expect_that(sol, equals(c(0.336737433210307, 0.663262566789693)))
})
