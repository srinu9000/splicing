
getExpr <- function(trueExpr, bioVar, noReads) {
  expr <- sapply(1:length(trueExpr), function(j) {
    rnorm(1, mean=trueExpr[j], sd=sqrt(bioVar[j]))
  })
  expr[expr < 0] <- 0
  expr/sum(expr)
}

#########################################################################

## Dominant isoform test.
## Is one isoform dominantly more abundant than the rest?
## I.e. can we say with confidence that a given isoform is
## expressed more highly than all other together?
## We do a one-sided t-test, against 1/2.

powDominant <- function(gene,           # gene structure
                        trueExpr,       # true expression levels (sums to 1)
                        domIso=1,       # the isoform to test
                        sampleSize=3,   # number of biological replicates
                        bioVar=0,       # biological variation, can be vector
                        readLength=30,  # read length
                        ## number of reads, defaults to 10 reads per base
                        noReads=sum(isoLength(gene)[[1]]*trueExpr)/30*10,
                        ## number of replicates for the simulation
                        noRep=100,
                        ## print dots to the console for each replicate
                        verbose=TRUE
                        ) {
  
  ## TODO: check arguments

  bioVar <- rep(bioVar, length.out=length(trueExpr))

  pvals <- sapply(1:noRep, function(i1) {
    reads <- lapply(1:sampleSize, function(i2) {
      expr <- getExpr(trueExpr, bioVar, noReads)
      simulateReads(gene, expression=expr, noReads=noReads,
                    readLength=readLength)
    })
    sols <- lapply(reads, function(r) { solveIso(gene, reads=r)$expression })
    isoToTest <- sapply(sols, "[", domIso)
    t.test(isoToTest, alternative="greater", mu=1/2)$p.value
  })
  
  pvals
}


###########################################################################

## Abundant isoform test
## Is one isoform significantly more abundant than another one?
## This is one-sided paired t-test with unequal variances

powAbund <- function(gene,              # gene structure
                     trueExpr,          # true expression levels (sums to 1)
                     majIso,            # the isoform to test, high expr.
                     minIso,            # the isoform to test, low expression
                     sampleSize=3,
                     bioVar=0,
                     readLength=30,
                     noReads=sum(isoLength(gene)[[1]]*trueExpr)/30*10,
                     ## number of replicates for the simulation
                     noRep=100,
                     ## print dots to the console for each replicate
                     verbose=TRUE
                     ) {

  ## TODO: check arguments

  bioVar <- rep(bioVar, length.out=length(trueExpr))

  pvals <- sapply(1:noRep, function(i1) {
    reads <- lapply(1:sampleSize, function(i2) {
      expr <- getExpr(trueExpr, bioVar, noReads)
      simulateReads(gene, expression=expr, noReads=noReads,
                    readLength=readLength)
    })
    sols <- lapply(reads, function(r) { solveIso(gene, reads=r)$expression })
    isosToTest <- sapply(sols, "[", c(majIso, minIso))
    t.test(isosToTest[1,], isosToTest[2,], alternative="greater",
           paired=TRUE)$p.value
  })
  
  pvals
}
  
###########################################################################

## Compare two samples

powCompSamp <- function(gene,
                        trueExpr1,
                        trueExpr2,
                        compIso,
                        sampleSize=3,
                        bioVar1=0,
                        bioVar2=0,
                        readLength=30,
                        noReads1=sum(isoLength(gene)[[1]]*trueExpr1)/30*10,
                        noReads2=sum(isoLength(gene)[[1]]*trueExpr2)/30*10,
                        noRep=100,
                        ## print dots to the console for each replicate
                        verbose=TRUE
                        ) {
  
  ## TODO: check arguments

  bioVar1 <- rep(bioVar1, length.out=length(trueExpr1))
  bioVar2 <- rep(bioVar2, length.out=length(trueExpr2))
  
  pvals <- sapply(1:noRep, function(i1) {
    reads1 <- lapply(1:sampleSize, function(i2) {
      expr <- getExpr(trueExpr1, bioVar1, noReads1)
      simulateReads(gene, expression=expr, noReads=noReads1,
                    readLength=readLength)
    })
    reads2 <- lapply(1:sampleSize, function(i2) {
      expr <- getExpr(trueExpr2, bioVar2, noReads2)
      simulateReads(gene, expression=expr, noReads=noReads2,
                    readLength=readLength)
    })
    
    sols1 <- lapply(reads1, function(r)
                    { solveIso(gene, reads=r)$expression })
    sols2 <- lapply(reads2, function(r)
                    { solveIso(gene, reads=r)$expression })

    iso1ToTest <- sapply(sols1, "[", compIso)
    iso2ToTest <- sapply(sols2, "[", compIso)
    
    t.test(iso1ToTest, iso2ToTest, alternative="greater",
           paired=TRUE)$p.value
  })
  
  pvals
}
  
###########################################################################

## Expression profile is different or not, in two samples
## We do a simple chi-square test

powProfile <- function(gene,
                       trueExpr1,
                       trueExpr2,
                       bioVar1=0,
                       bioVar2=0,
                       readLength=30,
                       noReads=sum(isoLength(gene)[[1]]*trueExpr)/30*10,
                       noRep=100,
                       ## print dots to the console for each replicate
                       verbose=TRUE
                       ) {

  ## TODO: check arguments

  bioVar1 <- rep(bioVar1, length.out=length(trueExpr1))
  bioVar2 <- rep(bioVar2, length.out=length(trueExpr2))
  
  pvals <- sapply(1:noRep, function(i1) {
    expr1 <- getExpr(trueExpr1, bioVar1, noReads)
    reads1 <- simulateReads(gene, expression=expr1, noReads=noReads,
                            readLength=readLength)
    expr2 <- getExpr(trueExpr2, bioVar2, noReads)
    reads2 <- simulateReads(gene, expression=expr2, noReads=noReads,
                            readLength=readLength)
    
    sol1 <- solveIso(gene, reads=reads1)$expression * noReads
    sol2 <- solveIso(gene, reads=reads2)$expression * noReads
    
    chisq.test(rbind(sol1, sol2))$p.value
  })
  
  pvals
}

###########################################################################

## Is an isoform expressed?
## One-sample t-test, comparing to zero

powPresent <- function(gene,
                       trueExpr,
                       checkIso=1,
                       sampleSize=3,
                       bioVar=0,
                       readLength=30,
                       noReads=sum(isoLength(gene)[[1]]*trueExpr)/30*10,
                       noRep=100,
                       ## print dots to the console for each replicate
                       verbose=TRUE
                       ) {

  ## TODO: check arguments

  bioVar <- rep(bioVar, length.out=length(trueExpr))

  pvals <- sapply(1:noRep, function(i1) {
    reads <- lapply(1:sampleSize, function(i2) {
      expr <- getExpr(trueExpr, bioVar, noReads)
      simulateReads(gene, expression=expr, noReads=noReads,
                    readLength=readLength)
    })
    sols <- lapply(reads, function(r) { solveIso(gene, reads=r)$expression })
    isoToTest <- sapply(sols, "[", checkIso)
    t.test(isoToTest, alternative="greater", mu=0)$p.value
  })
  
  pvals
}

