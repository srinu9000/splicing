\name{Splicing complexity}
\alias{crComplexity}
\title{Complexity of alternative splicing}
\description{
  TODO
}
\usage{
crComplexity (geneStructure, gene = 1, readLength, overHang = 1L, 
    expr, paired = FALSE, fast = FALSE, fragmentProb = NULL, 
    fragmentStart = 0L, normalMean = NA, normalVar = NA, numDevs = 4) 
}
\arguments{
  \item{geneStructure}{TODO}
  \item{gene}{TODO}
  \item{readLength}{TODO}
  \item{overHang}{TODO}
  \item{paired}{TODO}
  \item{fast}{TODO}  
  \item{fragmentProb}{TODO}
  \item{fragmentStart}{TODO}
  \item{normalMean}{TODO}
  \item{normalVar}{TODO}
  \item{numDevs}{TODO}
}
\details{
  TODO
}
\value{
  TODO
}
% \references{TODO}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \seealso{TODO}
\examples{
## TODO
my_gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                      list(c(1,2,3), c(1,3), c(1,2)))
crComplexity(my_gene, readLength=33, expr=c(1/3, 1/3, 1/3))

my_genes <- lapply(1:100, function(l) {
  createGene(list(c(1,100), c(201,200+l), c(401,500)),
             list(c(1,2,3), c(1,3)))
})
sapply(my_genes, crComplexity, readLength=33, expr=c(1/2, 1/2))
}
\keyword{misc}
