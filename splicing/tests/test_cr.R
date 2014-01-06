
context("CR complexity")

test_that("CR complexity works for exp vectors", {

  library(splicing)  
  gene3 <- createGene(list(c(1,100), c(201,700)-50, c(801,900)-100,
                         c(801,850)-100),
                      list(c(1,2,3), c(1,3), c(1,4)), id="g-3")
  
  cr <- crComplexity(gene3, readLength=50, expr=c(1,1,1)/3)
  cr0 <- structure(c(0.382706988947236, 4.53615969307819, 5.1066457027542),
                   .Dim = c(3L, 1L))
  test_that(cr, equals(cr0))
})

test_that("CR complexity works for exp matrices", {

  library(splicing)
  gene3 <- createGene(list(c(1,100), c(201,700)-50, c(801,900)-100,
                         c(801,850)-100),
                      list(c(1,2,3), c(1,3), c(1,4)), id="g-3")

  expr <- cbind(e1=c(1,1,1)/3, e2=c(1,2,3)/6)
  cr <- crComplexity(gene3, readLength=50, expr=expr)
  cr0 <- structure(c(0.382706988947236, 4.53615969307819, 5.1066457027542,
                     0.104846935140322, 2.39054500697885, 2.5700930134488),
                   .Dim = c(3L, 2L), .Dimnames = list(NULL, c("e1", "e2")))
  test_that(cr, equals(cr0))
})

test_that("CR complexity works when exp is not given", {

  library(splicing)
  set.seed(42)
  gene3 <- createGene(list(c(1,100), c(201,700)-50, c(801,900)-100,
                         c(801,850)-100),
                      list(c(1,2,3), c(1,3), c(1,4)), id="g-3")

  cr <- crComplexity(gene3, readLength=50)
  crm <- rowMeans(cr)
  crm0 <- c(0.503585690069611, 5.32086007756482, 5.87112677956185)
  test_that(crm, equals(crm0))
})
