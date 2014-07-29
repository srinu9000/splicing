
context("Power calculation for pooled samples")

test_that("Pooled power works for a single isoform", {

  set.seed(42)
  library(splicing)
  library(digest)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_pooled(gene, isoform = 1, coverage = 2)
  pow2 <- power_pooled(gene, isoform = 1, coverage = 10)

  expect_that(digest(pow1), equals("ea4a655a7b1ab489bbca8e762fa313ca"))
  expect_that(digest(pow2), equals("cc9cc06d672643435ac49681dbf9bd6e"))
  
})

test_that("Pooled power works for all isoforms", {

  set.seed(42)
  library(splicing)
  library(digest)
  
  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  pow1 <- power_pooled_all_iso(gene, coverage = 2)
  pow2 <- power_pooled_all_iso(gene, coverage = 10)

  expect_that(digest(pow1), equals("81fde082b0ce2838d33ac75ebc57be4f"))
  expect_that(digest(pow2), equals("7a9292147c15801b2bc62a22a6f8139f"))
  
})
