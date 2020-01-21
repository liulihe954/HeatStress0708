BiocManager::install("bladderbatch")
library(bladderbatch)
data(bladderdata)

library(limma)
BiocManager::install(c("sva","pamr"))
library(sva);library(pamr)
BiocManager::install(c("recount"))
library("recount", quietly = T);find.package("recount")
library("sva")
library("recount", quietly = T)

# library("WGCNA", quietly = T)



q_normalize <- function(dat){
  n = nrow(dat)
  p = ncol(dat)
  rank.dat =  dat # matrix for ranking
  for (i in 1:p){
    rank.dat[,i] = rank(dat[,i])
  }
  U = rank.dat/(n+1)
  qnorm(U)
}
Correct_pca = function(rse_raw,method){
  rse_raw <- t(rse_raw)# transpose data so that rows are samples and columns are gene expression measurements
  mod=matrix(1,nrow=dim(rse_raw)[1],ncol=1)
  colnames(mod)="Intercept"
  ## num.sv requires data matrix with features(genes) in the rows and samples in the column
  nsv=num.sv(t(rse_raw), mod, method = method)
  print(paste("Number of PCs estimated to be removed:", nsv))
  ## PC residualization of gene expression data using sva_network. Rows correspond to samples, Columns correspond to features
  exprs_corrected = sva_network(rse_raw, nsv)
  ## Quantile normalize the corrected gene expression measurements such that each expression of every gene follows a gaussian distribution
  exprs_corrected_norm <- q_normalize(exprs_corrected)
  return(list(exprs_corrected_norm = exprs_corrected_norm))
}


test_correction = Correct_pca(rse_raw,"leek")

test_correction$exprs_corrected_norm


originwd = getwd()
setwd("/Users/liulihe95/Desktop/HeatStress0708/networks_correction-master/publication_rmd/tutorial_vignette")
load("liver_subset_gtex.Rdata")
setwd(originwd);getwd()


dim(rse_raw);head(rse_raw)
dim(rse_raw)
rse_raw <- t(rse_raw) # transpose data so that rows are samples and columns are gene expression measurements
mod=matrix(1,nrow=dim(rse_raw)[1],ncol=1)
colnames(mod)="Intercept"
nsv=num.sv(t(rse_raw), mod, method = "leek") ## num.sv requires data matrix with features(genes) in the rows and samples in the column
print(paste("Number of PCs estimated to be removed:", nsv))
## PC residualization of gene expression data using sva_network. Rows correspond to samples, Columns correspond to features
exprs_corrected = sva_network(rse_raw, nsv)
## Quantile normalize the corrected gene expression measurements such that each expression of every gene follows a gaussian distribution
exprs_corrected_norm <- q_normalize(exprs_corrected)
dim(exprs_corrected_norm)

####
originwd = getwd()
setwd("/Users/liulihe95/Desktop/HeatStress0708/networks_correction-master/publication_rmd/tutorial_vignette")
load("rse_gene_liver.rdata")
setwd(originwd);getwd()


dim(rse_gene)
head(rse_gene)


rse <- scale_counts(rse_gene)

rse_raw <- log2(rse@assays$data$counts+2)



head(rse_raw)

genes_var <- apply(rse_raw, 1, var)
select_genes <- names(genes_var)[order(genes_var, decreasing = T)][1:1000]
rse_raw <- rse_raw[select_genes,]
save(rse_raw, file = "liver_subset_gtex.Rdata")




