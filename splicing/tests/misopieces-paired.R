
library(splicing)

options(width=60)

## REASSIGN_SAMPLES_PAIRED

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

rl <- 35
reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=rl, paired=TRUE,
                       normalMean=rl*2, normalVar=1, numDevs=0)

matches <- matchIso(gene, reads=reads, normalMean=rl*2,
                    normalVar=1, numDevs=4)
match.order <- order(apply(matches[[1]], 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
noiso <- 3L

set.seed(42)
r1 <- .Call("R_splicing_reassign_samples_paired", matches[[1]],
            match.order, psi, noiso, 1L, 1L, PACKAGE="splicing")

sec <- seq(1, noPairs(reads)*2, by=2)
calcCigar <- function(r) {
  fl <- r$position[sec+1] - r$position[sec]
  gr <- grepl("N", r$cigar[sec])
  sep <- ifelse(fl==rl | gr, "", paste(fl-rl, sep="", "N"))
  paste(sep="", r$cigar[sec], sep, r$cigar[sec+1])
}

reads2 <- list(chrname=reads$chrname, chrlen=reads$chrlen,
               chr=reads$chr[sec], qname=reads$qname[sec],
               cigar=calcCigar(reads),
               position=reads$position[sec], flag=rep(0L, noPairs(reads)),
               pairpos=rep(0L, noPairs(reads)), noPairs=0L,
               noSingles=noPairs(reads), paired=FALSE, mapq=reads$mapq[sec],
               rnext=rep(0L, noPairs(reads)), tlen=rep(0L, noPairs(reads)),
               seq=rep("*", noPairs(reads)), qual=rep("*", noPairs(reads)),
               mypair=rep(-1L, noPairs(reads)),
               attributes=reads$attributes[sec],
               sampleProb=reads$sampleProb)
class(reads2) <- "splicingSAM"

matches2 <- matchIso(gene, reads=reads2)

all((matches2 == (matches[[1]] != 0) + 0))

set.seed(42)
r2 <- .Call("R_splicing_reassign_samples", matches2,
            match.order, psi, noiso, 1L, 1L,
            PACKAGE="splicing")

all(r1 ==r2)

## SCORE_ISO_PAIRED

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

rl <- 35
reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=rl, paired=TRUE,
                       normalMean=100, normalVar=100, numDevs=4)

matches <- matchIso(gene, reads=reads, normalMean=100,
                    normalVar=100, numDevs=4)
match.order <- order(apply(matches[[1]], 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
noiso <- 3L
assignment <- .Call("R_splicing_reassign_samples_paired",
                    matches[[1]], match.order, psi, noiso, 1L, 1L,
                    PACKAGE="splicing")

hyper <- rep(1, noIso(gene))
isolen <- isoLength(gene)[[1]]
dist <- .Call("R_splicing_normal_fragment", 100, 100, 4, 2L*35L,
              PACKAGE="splicing")
lp <- sum(dist$fragmentStart - seq_along(dist$fragmentProb) + 2)
assscores <- log(isolen * length(dist$fragmentProb) - lp)

.Call("R_splicing_score_iso_paired", psi, 3L, assignment, isolen,
      assscores, PACKAGE="splicing")

psi2 <- cbind(c(5/10, 3/10, 2/10))
assignment <- .Call("R_splicing_reassign_samples_paired",
                    matches[[1]], match.order, psi2, noiso, 1L, 1L,
                    PACKAGE="splicing")

.Call("R_splicing_score_iso_paired", psi, 3L, assignment, isolen,
      assscores, PACKAGE="splicing")

## SCORE_JOINT_PAIRED

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

rl <- 35
reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=rl, paired=TRUE,
                       normalMean=100, normalVar=100, numDevs=4)

matches <- matchIso(gene, reads=reads, normalMean=100,
                    normalVar=100, numDevs=4)
match.order <- order(apply(matches[[1]], 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
noiso <- 3L
assignment <- .Call("R_splicing_reassign_samples_paired",
                    matches[[1]], match.order, psi, noiso, 1L, 1L,
                    PACKAGE="splicing")

hyper <- rep(1, noIso(gene))
isolen <- isoLength(gene)[[1]]
dist <- .Call("R_splicing_normal_fragment", 100, 100, 4, 2L*35L,
              PACKAGE="splicing")
noiso <- 3L
matches <- matchIso(gene, reads=reads, normalMean=100,
                    normalVar=100, numDevs=4)
fragmentLength <- matches[[2]]

isoscores <- matrix(0, nr=length(dist$fragmentProb), nc=noiso)
assscores <- numeric(noiso)
for (j in 1:length(dist$fragmentProb)) {
  logprob <- dist$fragmentProb[j]
  for (i in 1:noiso) {
    lp <- isolen[i] - dist$fragmentStart - j + 1
    isoscores[j,i] <- -log(lp) + logprob
    assscores[i] <- assscores[i] + lp
  }
}
assscores <- log(assscores)

.Call("R_splicing_score_joint_paired", assignment,
      as.integer(noPairs(reads)), 1L, psi, hyper,
      isolen, isoscores, assscores, fragmentLength,
      dist$fragmentStart, PACKAGE="splicing")

psi2 <- cbind(c(5/10, 3/10, 2/10))
assignment <- .Call("R_splicing_reassign_samples_paired",
                    matches[[1]], match.order, psi2, noiso, 1L, 1L,
                    PACKAGE="splicing")

.Call("R_splicing_score_joint_paired", assignment,
      as.integer(noPairs(reads)), 1L, psi, hyper,
      isolen, isoscores, assscores, fragmentLength,
      dist$fragmentStart, PACKAGE="splicing")

## METROPOLIS_HASTINGS_RATIO_PAIRED

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

rl <- 35
reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=rl, paired=TRUE,
                       normalMean=100, normalVar=100, numDevs=4)

matches <- matchIso(gene, reads=reads, normalMean=100,
                    normalVar=100, numDevs=4)
match.order <- order(apply(matches[[1]], 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
noiso <- 3L
assignment <- .Call("R_splicing_reassign_samples_paired",
                    matches[[1]], match.order, psi, noiso, 1L, 1L,
                    PACKAGE="splicing")
alpha <- cbind(c(1/3, 2/3))

psinew <- cbind(c(5/10, 3/10, 2/10))
alphanew <- cbind(c(1/2, 1/2))
sigma <- 0.05
hyperp <- rep(1, noiso)

isolen <- isoLength(gene)[[1]]
dist <- .Call("R_splicing_normal_fragment", 100, 100, 4, 2L*35L,
              PACKAGE="splicing")
fragmentLength <- matches[[2]]

isoscores <- matrix(0, nr=length(dist$fragmentProb), nc=noiso)
assscores <- numeric(noiso)
for (j in 1:length(dist$fragmentProb)) {
  logprob <- dist$fragmentProb[j]
  for (i in 1:noiso) {
    lp <- isolen[i] - dist$fragmentStart - j + 1
    isoscores[j,i] <- -log(lp) + logprob
    assscores[i] <- assscores[i] + lp
  }
}
assscores <- log(assscores)

.Call("R_splicing_metropolis_hastings_ratio_paired", assignment,
      as.integer(noPairs(reads)), 1L,
      psinew, alphanew, psi, alpha, sigma, noiso, isolen, hyperp,
      isoscores, assscores, matches[[2]], dist$fragmentStart, 1L,
      PACKAGE="splicing")

.Call("R_splicing_metropolis_hastings_ratio_paired", assignment,
      as.integer(noPairs(reads)), 1L,
      psi, alpha, psinew, alphanew, sigma, noiso, isolen, hyperp,
      isoscores, assscores, matches[[2]], dist$fragmentStart, 1L,
      PACKAGE="splicing")

