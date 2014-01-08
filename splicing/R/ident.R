
assignmentMatrix <- function(geneStructure, gene=1L, readLength,
                             overHang=1L, bias=0, paired=FALSE, fast=FALSE,
                             fragmentProb=NULL, fragmentStart=0L,
                             normalMean=NA, normalVar=NA, numDevs=4) {

  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
  
  if (!paired) { 
    res <- .Call("R_splicing_assignment_matrix_bias", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 as.integer(overHang), as.integer(bias),
                 PACKAGE="splicing")
  } else {
    res <- .Call("R_splicing_paired_assignment_matrix", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 as.integer(overHang), as.integer(bias),
                 as.logical(fast), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 PACKAGE="splicing")
  }
  
  colnames(res) <- apply(res, 2, function(x) {
    paste(ifelse(x==0, '0', '1'), collapse="")
  })
  res
}

geneComplexity <- function(geneStructure, gene=1, readLength, overHang=1L,
                           type=c("relative", "absolute"),
                           norm=c("2","1","inf"), paired=FALSE, fast=FALSE,
                           fragmentProb=NULL, fragmentStart=0L, normalMean=NA,
                           normalVar=NA, numDevs=4) {

  type <- switch(match.arg(type), "relative"=0, "absolute"=1)
  norm <- switch(match.arg(norm), "2"=0, "1"=1, "inf"=2)
  if (!is.null(fragmentProb)) { fragmentProb <- as.numeric(fragmentProb) }

  .Call("R_splicing_gene_complexity", geneStructure, as.integer(gene),
        as.integer(readLength), as.integer(overHang), as.integer(type),
        as.integer(norm), as.logical(paired), as.logical(fast), fragmentProb,
        as.integer(fragmentStart), as.double(normalMean),
        as.double(normalVar), as.double(numDevs),
        PACKAGE="splicing")
} 

isoComplexity <- function(geneStructure, gene=1, readLength, overHang=1L,
                          expr, noReads=1,
                          paired=FALSE, fast=FALSE, fragmentProb=NULL,
                          fragmentStart=0L, normalMean=NA, normalVar=NA,
                          numDevs=4) {

  require(MASS)

  if (missing(expr)) {
    noi <- noIso(geneStructure)[gene]
    expr <- rep(1/noi, noi)
  }

  mat <- assignmentMatrix(geneStructure, gene=gene, readLength=readLength,
                          overHang=overHang, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  mat <- t(mat / rowSums(mat))

  il <- isoLength(geneStructure)[[gene]]
  nf <- sum(il * expr)
  e2 <- il * expr / nf

  cc <- as.vector(mat %*% e2)
  V <- -outer(cc, cc)
  diag(V) <- cc * (1-cc)
  V <- V * noReads

  matInv <- ginv(mat)
  V2 <- matInv %*% V %*% t(matInv) / noReads^2
  nf2 <- sum(e2 / il)
  V3 <- V2 / nf2^2
  V3 <- V3 / il
  V3 <- t(t(V3)/il)
  diag(V3)
}

crComplexity <- function(geneStructure, gene=1, readLength, overHang=1L,
                         expr, paired=FALSE, fast=FALSE,
                         fragmentProb=NULL, fragmentStart=0L, normalMean=NA,
                         normalVar=NA, numDevs=4) {

  require(MASS)
  
  rdirichlet <- function (n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
  }

  noi <- noIso(geneStructure)[gene]
  if (missing(expr)) {
    expr <- t(rdirichlet(100, rep(1, noi)))
  } else {
    if (is.matrix(expr)) {
      expr <- cbind(expr)
    } else {
      expr <- cbind(expr)
      colnames(expr) <- NULL
    }
  }
  
  mycr <- function(expr) {
  
    mat <- assignmentMatrix(geneStructure, gene=gene, readLength=readLength,
                            overHang=overHang, expr, paired=paired, fast=fast,
                            fragmentProb=fragmentProb,
                            fragmentStart=fragmentStart,
                            normalMean=normalMean, normalVar=normalVar,
                            numDevs=numDevs)
    
    elen <- rowSums(mat)
    elen2 <- isoLength(geneStructure)[[gene]] - readLength + 1
    if (any(elen != elen2)) { stop("Probably wrong gene structure") }
    mat <- t(mat / elen)
    noc <- nrow(mat)
    
    X <- mat %*% (elen * expr)
    Y <- sum(elen * expr)
    
    I <- matrix(0, noi-1, noi-1)
    for (k1 in 1:(noi-1)) {
      for (k2 in k1:(noi-1)) {
        C <- 0
        for (i in 1:noc) {
          if (X[i] > 1e-15) {
            t1 <- (elen[k1]-elen[noi]) * (elen[k2]-elen[noi]) / (Y^2)
            t2 <- (mat[i,k1]*elen[k1] - mat[i,noi]*elen[noi]) *
              (mat[i,k2]*elen[k2] - mat[i,noi]*elen[noi]) / (X[i]^2)
            C <- C + (X[i]/Y) * (t1 - t2)
          }
        }
        I[k1,k2] <- I[k2,k1] <- -C
      }
    }
    
    I[I==Inf] <- 1e100
    cr <- ginv(I)
    c(diag(cr), sum(cr))
  }

  cr <- apply(expr, 2, mycr)
  colnames(cr) <- colnames(expr)
  cr
}
