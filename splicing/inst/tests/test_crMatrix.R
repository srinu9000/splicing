
context("Fisher information matrix")

test_that("crMatrix function works", {

  library(splicing)

  gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                    list(c(1,2,3), c(1,3)))

  M1 <- crMatrix(gene, readLength=10, expr=c(.5,.5))
  M2 <- crMatrix(gene, readLength=50, expr=c(.5,.5))
  M3 <- crMatrix(gene, readLength=100, expr=c(.5,.5))
  ## TODO: M4 <- crMatrix(gene, readLength=201, expr=c(.5,.5))

  d.M1 <- structure(c(1.23903490058879, -1.23903490058879,
                     -1.23903490058879, 1.23903490058879),
                   .Dim = c(2L, 2L))
  d.M2 <- structure(c(0.580507500431059, -0.580507500431059,
                     -0.580507500431059, 0.580507500431059),
                   .Dim = c(2L, 2L))
  d.M3 <- structure(c(0.285026751337567, -0.285026751337567,
                     -0.285026751337567, 0.285026751337567),
                   .Dim = c(2L, 2L))
  
  expect_that(M1, equals(d.M1))
  expect_that(M2, equals(d.M2))
  expect_that(M3, equals(d.M3))

  gene2 <- createGene(list(c(1,100), c(201,300), c(401,500)),
                     list(c(1,2,3), c(1,3), c(1,2)))

  M1.2 <- crMatrix(gene2, readLength=10, expr=c(1,1,1)/3)
  M2.2 <- crMatrix(gene2, readLength=50, expr=c(1,1,1)/3)
  M3.2 <- crMatrix(gene2, readLength=100, expr=c(1,1,1)/3)

  d.M1.2 <- structure(c(3.77091257108143, -1.64188059315029, -2.12903197793114, 
                       -1.64188059315029, 1.43046262934082, 0.211417963809467,
                       -2.12903197793114, 0.211417963809467, 1.91761401412168),
                     .Dim = c(3L, 3L))
  d.M2.2 <- structure(c(1.20022919695949, -0.245830076485678, -0.95439912047381, 
                       -0.245830076485678, 0.543880412809873, -0.298050336324195,
                       -0.95439912047381, -0.298050336324195, 1.252449456798),
                     .Dim = c(3L, 3L))
  d.M3.2 <- structure(c(0.599559728742827, -0.00200521648408972, -0.597554512258738, 
                       -0.00200521648408972, 0.299523614622206, -0.297518398138116, 
                       -0.597554512258738, -0.297518398138116, 0.895072910396854),
                     .Dim = c(3L, 3L))

  
  expect_that(M1.2, equals(d.M1.2))
  expect_that(M2.2, equals(d.M2.2))
  expect_that(M3.2, equals(d.M3.2))

})
