library(MeSH.db)
library(MeSH.Bta.eg.db)
library(meshr)

showMethods("category")
showMethods("database")
data(geneid.cummeRbund)
names(geneid.cummeRbund)

data(geneid.cummeRbund)
names(geneid.cummeRbund)
## This data is also available by following scripts.
if(interactive()){
  library(cummeRbund)
  library(org.Hs.eg.db)
  cuff <- readCufflinks(dir = system.file("extdata", package = "cummeRbund"))
  gene.symbols <- annotation(genes(cuff))[,4]
  mySigGeneIds <- getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='genes')
  mySigGenes <- getGenes(cuff,mySigGeneIds)
  sig.gene.symbols <- annotation(mySigGenes)[,4]
  gene.symbols <- gene.symbols[!is.na(gene.symbols)]
  sig.gene.symbols <- sig.gene.symbols[!is.na(sig.gene.symbols)]
  geneid.cummeRbund <- select(org.Hs.eg.db, keys=gene.symbols, keytype="SYMBOL", columns="ENTREZID")
  sig.geneid.cummeRbund <- select(org.Hs.eg.db, keys=sig.gene.symbols, keytype="SYMBOL", columns="ENTREZID")
  na.index1 <- which(is.na(geneid.cummeRbund[,2]))
  for (i in na.index1){
    s <- unlist(strsplit(as.character(geneid.cummeRbund[i,][1]), ","))[1]
    sym <- get(s, org.Hs.egALIAS2EG)[1]
    geneid.cummeRbund[i,2] <- as.integer(sym)
  }
  na.index2 <- which(is.na(sig.geneid.cummeRbund[,2]))
  for (i in na.index2){
    s <- unlist(strsplit(as.character(sig.geneid.cummeRbund[i,][1]), ","))[1]
    sym <- get(s, org.Hs.egALIAS2EG)[1]
    sig.geneid.cummeRbund[i,2] <- as.integer(sym)
  }
  geneid.cummeRbund <- geneid.cummeRbund[!duplicated(geneid.cummeRbund[,2]), ]
  sig.geneid.cummeRbund <- sig.geneid.cummeRbund[!duplicated(sig.geneid.cummeRbund[,2]), ]
}

