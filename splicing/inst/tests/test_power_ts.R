
context("Power for trend in time series")


test_that("Power for trend in time series, one isoform", {

  set.seed(42)
  library(splicing)
  library(digest)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_ts(gene, isoform = 1, coverage = 2)
  pow2 <- power_ts(gene, isoform = 1, coverage = 10)

  expect_that(digest(pow1), equals("18a0c881fcd0d77afbb50f855e0e8b9b"))
  expect_that(digest(pow2), equals("713df9b65f23e1f36d83d0a09eeffe30"))

})

test_that("Power for trend in time series, all isoforms", {

  set.seed(42)
  library(splicing)
  library(digest)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_ts_all_iso(gene, coverage = 2)
  pow2 <- power_ts_all_iso(gene, coverage = 10)

  expect_that(digest(pow1), equals("b0e993fc7aaef5e5716a25e1e8e80bc7"))
  expect_that(digest(pow2), equals("5fe1da44b2b50e1a79438d669d096bdb"))

})
