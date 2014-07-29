
context("Power for comparing two groups of samples")

test_that("Power for two groups works for a single isoform", {

  set.seed(42)
  library(splicing)
  library(digest)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_rep(gene, isoform = 1, coverage = 2)
  pow2 <- power_rep(gene, isoform = 1, coverage = 10)

  expect_that(digest(pow1), equals("fb7e548dec35d15aafbb19754a7d07ee"))
  expect_that(digest(pow2), equals("7821a5a5e44e370203b6f716d4e9221f"))

})

test_that("Power, two groups, all isforms", {

  set.seed(42)
  library(splicing)
  library(digest)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_rep_all_iso(gene, coverage = 2)
  pow2 <- power_rep_all_iso(gene, coverage = 10)

  expect_that(digest(pow1), equals("930265f4bbe1169483d746219d1e29c0"))
  expect_that(digest(pow2), equals("4b02e966b63afc8b8c064706cbf50b7a"))

})
