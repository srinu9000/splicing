\name{Power simulations}
\alias{powAbund}
\alias{powCompSamp}
\alias{powDominant}
\alias{powPresent}
\alias{powProfile}
\title{Power simulations for isoform quantitation}
\description{
  TODO
}
\usage{
powAbund(gene, trueExpr, majIso, minIso, sampleSize=3, bioVar=0,
         readLength=30, noReads=sum(isoLength(gene)[[1]] * trueExpr)/30 * 10,
         noRep=100, verbose=TRUE)
}
\arguments{
  \item{}{}
}
\details{
  TODO
}
\value{
  TODO
}
\references{TODO}
\author{Gabor Csardi \email{gcsardi@stat.harvard.edu}}
\seealso{TODO}
\examples{
## Two genes, a simple and one that is hard to estimate
g1 <- createGene(list(c(1L,1000L), c(901L,1900L)), list(1L,2L))
g2 <- createGene(list(c(1L,1000L), c(51L,1050L)), list(1L,2L))

dev.new(width=12, height=4)
layout(matrix(1:12, nrow=2))
plotIso(g1, mar=c(5,4,4,2)+.1)
plotIso(g2, mar=c(5,4,4,2)+.1)

## Dominant isoforms
p11 <- powDominant(g1, trueExpr=c(2/3, 1/3))
p21 <- powDominant(g2, trueExpr=c(2/3, 1/3))

hist(p11, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Simple gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")
hist(p21, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Complex gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")

## More abundant isoform
p12 <- powAbund(g1, trueExpr=c(2/3, 1/3), majIso=1, minIso=2)
p22 <- powAbund(g2, trueExpr=c(2/3, 1/3), majIso=1, minIso=2)

hist(p12, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Simple gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")
hist(p22, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Complex gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")

## Compare the same isoform in two samples
p13 <- powCompSamp(g1, c(2/3, 1/3), c(1/3, 2/3), compIso=1)
p23 <- powCompSamp(g2, c(2/3, 1/3), c(1/3, 2/3), compIso=1)

hist(p13, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Simple gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")
hist(p23, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Complex gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")

## Compare two profiles
p14 <- powProfile(g1, c(2/3, 1/3), c(1/3, 2/3), noReads=100)
p24 <- powProfile(g2, c(2/3, 1/3), c(1/3, 2/3), noReads=100)

hist(p14, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Simple gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")
hist(p24, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Complex gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")

## Is an isoform expressed?
p15 <- powPresent(g1, trueExpr=c(1/3, 2/3), noReads=100)
p25 <- powPresent(g2, trueExpr=c(1/3, 2/3), noReads=100)

hist(p15, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Simple gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")
hist(p25, breaks=seq(0,1,by=.05), ylim=c(0,100), main="Complex gene",
     col=c("#00ff0099", rep("#00000099", 19)), xlab="p-values")

}
\keyword{misc}
