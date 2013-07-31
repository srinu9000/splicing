
na0 <- function(x) { x[is.na(x)] <- 0; x }

droplast <- function(x) { x[-length(x)] }

normlogit <- function(x) {
  x[x < 0] <- 1
  x <- x/sum(x)
  log(droplast(x)) - log(x[length(x)])
}

logit_inv <- function(x) {
  ex <- exp(x)
  if (any(ex == Inf)) {
    c(ifelse(ex == Inf, 1/sum(ex==Inf), 0), 0)
  } else {
    y <- exp(x) / (sum(exp(x))+1)
    c(y, 1 - sum(y))
  }
}

mleIso <- function(geneStructure, gene=1, reads,
                   readLength=getReadLength(reads), method="BFGS", ...) {

  require(MASS)
  
  ass <- assignmentMatrix(geneStructure, gene=gene,
                          readLength=readLength)
  ass <- t(ass / rowSums(ass))
  mat <- matchIso(geneStructure, gene=gene, reads=reads)
  matvec <- (function(m) {
    mv <- table(apply((m != 0) + 0L, 2, paste, collapse=""))[rownames(ass)]
    mv[is.na(mv)] <- 0
    names(mv) <- rownames(ass)
    mv
  })(mat)
  matvec <- as.vector(matvec)
  elen <- isoLength(geneStructure)[[gene]] - readLength + 1
  
  f <- function(le) {
    e <- logit_inv(le)
    nc <- sum(elen * e)
    e2 <- elen * e / nc
    ni <- ass %*% e2
    ni[ni < 0] <- 0                     # numerical errors, result is -Inf 
    res <- sum(matvec * as.vector(log(ni)), na.rm=TRUE)
    if (res == -Inf) { res <- -1e+50 }
    res
  } 
  ores <- optim(na0(normlogit(ginv(ass) %*% matvec)), f,
                control=list(fnscale=-1), method=method, ...)
  logit_inv(ores$par)
}
