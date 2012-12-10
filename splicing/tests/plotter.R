
gene <- createGene(list(c(1,200), c(501,600), c(701,850)),
                   list("iso1"=c(1,2,3), "iso2"=c(1,3)), id="insil")
gene$ID <- sub("-isoform-", "-iso-", gene$ID)
reads1 <- simulateReads(gene, readLength=30, noReads=1000,
                        expression=c(0.7,0.3))
reads2 <- simulateReads(gene, readLength=30, noReads=1000,
                        expression=c(.2, .8))
reads3 <- simulateReads(gene, readLength=30, noReads=1000,
                        expression=c(.5, .5))
reads4 <- simulateReads(gene, readLength=30, noReads=1000,
                        expression=c(.01, .99))
reads <- list(reads1=reads1, reads2=reads2, reads3=reads3, reads4=reads4)

misoResult <- lapply(reads, MISO, geneStructure=gene, gene=1L)

plotReads(gene, reads, misoResult)

### 

plotReads(gene, reads[-1], misoResult[-1])
plotReads(gene, reads[-(1:2)], misoResult[-(1:2)])

##

library(splicing)

set.seed(42)

gene <- createGene(list(c(101,200), c(501,600), c(750,900), c(1001, 1400)),
                   list(c(1,2,3), c(1,3), c(1,2,4)))

## size <- plotIsoSize(gene)
## dev.new(width=size["width"], height=size["height"])
## plotIso(gene)

reads1 <- simulateReads(gene, readLength=30, noReads=400,
                        expression=c(0.7,0.1,0.1,0.1))
reads2 <- simulateReads(gene, readLength=30, noReads=400,
                        expression=c(.1, .6,.1,.1))
reads3 <- simulateReads(gene, readLength=30, noReads=400,
                        expression=c(.25, .25, .25, .25))
reads4 <- simulateReads(gene, readLength=30, noReads=400,
                        expression=c(.01, .9, .05, .04))
reads <- list(reads1=reads1, reads2=reads2, reads3=reads3, reads4=reads4)

misoResult <- lapply(reads, MISO, geneStructure=gene, gene=1L)

dev.new(width=5, height=3)
plotMISO(misoResult[[1]], type="area")

size <- plotReadsSize(gene, reads, misoResult)
dev.new(width=size["width"], height=size["height"])
plotReads(gene, reads, misoResult)

