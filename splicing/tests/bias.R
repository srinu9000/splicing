
library(splicing)

set.seed(42)
options(width=60)

##################################################
## A gene for which the bias cannot be estimated

gene <- createGene(list(c(1,100), c(1,300), c(201,300)),
                   list(c(1,3), 2))

assignmentMatrix(gene, gene=1L, readLength=30)
cbind("01"=c(0,200-71),
      "10"=c(100-71,0),
      "11"=c(71-0+171-100,71-0+271-200))

assignmentMatrix(gene, gene=1L, readLength=30, bias=1L)
cbind("01"=c(0,200^2-71^2),
      "10"=c(100^2-71^2,0),
      "11"=c(71^2-0^2+171^2-100^2,71^2-0^2+271^2-200^2))

assignmentMatrix(gene, gene=1L, readLength=30, bias=2L)
cbind("01"=c(0,200^3-71^3),
      "10"=c(100^3-71^3,0),
      "11"=c(71^3-0^3+171^3-100^3,71^3-0^3+271^3-200^3))

reads <- simulateReads(gene, readLength=30, noReads=2000,
                       expression=c(0.8,0.2))

s1 <- solveIso(gene, reads=reads)
tryres <- try(s2 <- solveIsoLinBias(gene, reads=reads), silent=TRUE)
inherits(tryres, "try-error")

#################################################
## Linear bias can be estimated, but not the
## quadratic bias

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(2,3), 1, 2, c(1,2,3)))

reads <- simulateReads(gene, readLength=30, noReads=2000,
                       expression=c(0.2, 0.3, 0.3, 0.1, 0.1, 0))

s1 <- solveIso(gene, reads=reads)
s2 <- solveIsoLinBias(gene, reads=reads)
tryres <- try(s3 <- solveIsoQuadBias(gene, reads=reads), silent=TRUE)
inherits(tryres, "try-error")

round(s1$expression, 4)
round(s2$expression, 4)
round(s2$a, 4)

#################################################
## Quadratic bias can be estimated here

gene <- createGene(list(c(1,100), c(201,300), c(401,500), c(601,700)),
                   list(c(1,2), c(1,3), c(2,3), c(2,4)))

reads <- simulateReads(gene, readLength=30, noReads=200,
                       expression=c(0.2, 0.3, 0.3, 0.1))

s1 <- solveIso(gene, reads=reads)
s2 <- solveIsoLinBias(gene, reads=reads)
s3 <- solveIsoQuadBias(gene, reads=reads)

round(s1$expression, 4)
round(s2$expression, 4)
round(s3$expression, 4)
round(s2$a, 4)
round(s3$a, 4)
round(s3$b, 4)

##################################################
## Linear bias, paired-end reads

library(splicing)

set.seed(1)

gene <- createGene(list(c(1,100), c(1,300), c(201,300)),
                   list(c(1,3), 2))

cbind("01"=c(0,200-71),
      "10"=c(100-71,0),
      "11"=c(71-0+171-100,71-0+271-200))

assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=TRUE,
                 normalMean=30, normalVar=1e-10)
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=FALSE,
                 normalMean=30, normalVar=1e-10)
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=TRUE,
                 normalMean=30, normalVar=60)
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=TRUE,
                 normalMean=30, normalVar=60, bias=1L)
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=FALSE,
                 normalMean=30, normalVar=60)

reads <- simulateReads(gene, readLength=30, noReads=1000,
                       expression=c(0.2, 0.8),
                       paired=TRUE, normalMean=30, normalVar=30)

cbind("01"=c(0,200^2-71^2),
      "10"=c(100^2-71^2,0),
      "11"=c(71^2-0^2+171^2-100^2,71^2-0^2+271^2-200^2))
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=TRUE,
                 normalMean=30, normalVar=30, bias=1L)
assignmentMatrix(gene, gene=1L, readLength=10, paired=TRUE, fast=TRUE,
                 normalMean=30, normalVar=30, bias=1L)

s1 <- solveIso(gene, reads=reads, normalMean=30, normalVar=30, fast=TRUE)
s2 <- solveIsoLinBias(gene, reads=reads, normalMean=30, normalVar=30,
                      fast=TRUE)
round(s1$expression,4)
round(s2$expression,4)
round(s2$a, 4)

sort(unique(abs(reads$position - reads$pairpos)))

s1$nomatch

bad <- unique(reads$position - reads$pairpos)

for (i in 1:100) {
  sam <- sample(1:(noReads(reads)/2), 20)
  reads2 <- selectReads(reads, c(2*(sam-1)+1, 2*sam))
  s2 <- solveIsoLinBias(gene, reads=reads2, normalMean=30, normalVar=60,
                        fast=TRUE)
  if (s2$a < 1/2) {
    bad <- setdiff(bad, unique(reads2$position - reads2$pairpos))
  }
}

###############################################################
## Linear bias, paired-end reads

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(2,3), 1, 2, c(1,2,3)))

reads <- simulateReads(gene, readLength=10, noReads=1000,
                       expression=c(0.2, 0.3, 0.3, 0.1, 0.1, 0),
                       paired=TRUE, normalMean=30, normalVar=1e-9)

s1 <- solveIso(gene, reads=reads, normalMean=30, normalVar=1e-9, fast=TRUE)
s2 <- solveIsoLinBias(gene, reads=reads, normalMean=30, normalVar=1e-9,
                      fast=TRUE)

round(s1$expression,4)
round(s2$expression,4)
round(s2$a, 4)
