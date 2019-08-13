##############################################
require(sva);require(WGCNA)
require(ppcor);require(dplyr)
require(edgeR);require(clusterProfiler)
require(ggplot2);require(magrittr)
require(biomaRt);require(gage);require(doParallel)
require(limma);require(recount);require(pamr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
CowsID_ht = c("6334","8514","8971","8867","8841","8966")# ID of heat group
CowsID_cl = c("8252","8832","8896","8983","8897","8862")# ID of cool group
# (all the IDs) 27570 genes; 12 cows; 3 time points 
networkData = read.csv("count_matrix.tsv",sep = "\t",header = T);rownames(networkData) = networkData[,1];networkData = networkData[,-1] 
dim(networkData)
## find index 
# day 14
names14_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos14 = grep(".14",names14_last3)
names14 = names(networkData)[pos14];pos14_ht = substr(names14,2,5) %in% CowsID_ht;pos14_cl = substr(names14,2,5) %in% CowsID_cl
# day 42
names42_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos42 = grep(".42",names42_last3)
names42 = names(networkData)[pos42];pos42_ht = substr(names42,2,5) %in% CowsID_ht;pos42_cl = substr(names42,2,5) %in% CowsID_cl
# day 84
names84_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos84 = grep(".84",names84_last3)
names84 = names(networkData)[pos84];pos84_ht = substr(names84,2,5) %in% CowsID_ht;pos84_cl = substr(names84,2,5) %in% CowsID_cl
#compose index
column_14_cl = pos14[pos14_cl];column_14_ht = pos14[pos14_ht]
column_42_cl = pos42[pos42_cl];column_42_ht = pos42[pos42_ht]
column_84_cl = pos84[pos84_cl];column_84_ht = pos84[pos84_ht]
# -4 in day14; -4 in day 42; -7 in day84 (think about which one to respect)
networkData14 = networkData[,c(column_14_cl,column_14_ht)]
networkData42 = networkData[,c(column_42_cl,column_42_ht)]
networkData84 = networkData[,c(column_84_cl,column_84_ht)]
dim(networkData14);dim(networkData42);dim(networkData84)



################################################################################
#### substitute dataset for multiple scripts ###################################
networkData14 = networkData84                ###################################
#### substitute dataset for multiple scripts ###################################
################################################################################


########################################################################################################################
# step 1 - filter out top 40% counts
## filter out top 40% counts # function established for future use
remove_filter = function(networkData,thres){
  ID_meanexpr1 = data.frame(names = rownames(networkData), mean = apply(networkData, MARGIN = 1,mean));
  ID_meanexpr2 = cbind(ID_meanexpr1,percent = ID_meanexpr1$mean/sum(ID_meanexpr1$mean))
  ID_meanexpr3 = ID_meanexpr2[order(ID_meanexpr2$mean,decreasing = T),]
  accumulative = numeric(nrow(ID_meanexpr3))
  for (i in c(1:nrow(ID_meanexpr3))){
    accum = sum(ID_meanexpr3$percent[1:i])
    accumulative[i] = accum
  }
  remove_pos = (length(which(accumulative <= thres))+1)
  remove_index = ID_meanexpr3$names[1:remove_pos]
  networkData_filter = networkData[!(rownames(networkData)%in%remove_index),]
  Results = list(remove_index=remove_index,networkData_filter = networkData_filter)
  return(Results)
}
networkData14_filter = remove_filter(networkData14,0.4)$networkData_filter
#networkData42_filter = remove_filter(networkData42,0.4)$networkData_filter
#networkData84_filter = remove_filter(networkData84,0.4)$networkData_filter
dim(networkData14_filter)#;dim(networkData42_filter);dim(networkData84_filter);

# step 2 - normalization (0s out and normalization)
#BiocManager::install("edgeR") 
# zeros out!
remove_index14 = which(rowSums(networkData14_filter) == 0);length(remove_index14)
#remove_index42 = which(rowSums(networkData42_filter) == 0);length(remove_index42)
#remove_index84 = which(rowSums(networkData84_filter) == 0);length(remove_index84)

# normalization
require(edgeR)
networkData14_nm1 = networkData14_filter[-remove_index14,];dim(networkData14_nm1)
networkData14_nmList = DGEList(counts = networkData14_nm1,group  = c(rep("CL",length(column_14_cl)),rep("HT",length(column_14_ht))))
networkData14_nm2 = calcNormFactors(networkData14_nmList)
networkData14_normalized_normfactors = networkData14_nm2$samples
networkData14_normalized = data.frame(networkData14_nm2$counts)
dim(networkData14_normalized)

#load("networkData14 norm and factors.RData")
# step 3 - log 2 trans
# log trans
networkData14_log2 = log2(networkData14_normalized+2)

# step 4 - filter out bottom 50% variation
# select most var
networkData14_log2$variance = apply(networkData14_log2,1,var)
networkData14_log2_50var = networkData14_log2[networkData14_log2$variance >= quantile(networkData14_log2$variance,c(.50)),] #50% most variable genes
networkData14_log2_50var$variance <- NULL
dim(networkData14_log2_50var)

# step 5 - pca correction 
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
  return(list(exprs_corrected_norm = t(data.frame(exprs_corrected_norm))))
}
networkData14_correction = Correct_pca(networkData14_log2_50var,"leek")
networkData14_final = data.frame(networkData14_correction$exprs_corrected_norm); 
names(networkData14_final) = names(networkData14_log2_50var)
# step 6 - select into two groups
datExpr14_cl = t(networkData14_final[,colnames(networkData14_final) %in% names(networkData)[column_84_cl] ])
datExpr14_ht = t(networkData14_final[,colnames(networkData14_final) %in% names(networkData)[column_84_ht] ])
datExpr14_cl = data.frame(datExpr14_cl);datExpr14_ht = data.frame(datExpr14_ht)
dim(datExpr14_cl);dim(datExpr14_ht)

##########
print("check na results bicor")load("CoolHeatDay14_modulePreservation bicor.RData")
load("CoolHeatday14 bicor.RData")
ref=1; test = 2
### print results - short version
statsObs14_b = cbind(mp14_b$quality$observed[[ref]][[test]][,-1], mp14_b$preservation$observed[[ref]][[test]][,-1])
statsZ14_b = cbind(mp14_b$quality$Z[[ref]][[test]][, -1], mp14_b$preservation$Z[[ref]][[test]][,-1]);
# Compare preservation to quality:
print(cbind(statsObs14_b[,c("medianRank.pres", "medianRank.qual")],
            signif(statsZ14_b[,c("Zsummary.pres", "Zsummary.qual")],2)))
### print results - full
Obs.PreservationStats14_b= mp14_b$preservation$observed[[ref]][[test]]
Z.PreservationStats14_b=mp14_b$preservation$Z[[ref]][[test]]

modColors14_b = rownames(Obs.PreservationStats14_b)
moduleSize14_b = Obs.PreservationStats14_b$moduleSize
# we will omit the grey module (background genes) and the gold module (random sample of genes)
selectModules = !(modColors14_b %in% c("grey", "gold"))
# Text labels for points
point.label14_b = modColors14_b[selectModules]
#Composite preservation statistics
medianRank14_b=Obs.PreservationStats14_b$medianRank.pres
Zsummary14_b=Z.PreservationStats14_b$Zsummary.pres
# length(Zsummary14_b)

#===========================================================================================
#                             11. Gene Ontology enrichment                                ##
#===========================================================================================
#database2 = useMart("ensembl")
#genome2 = useDataset("btaurus_gene_ensembl", mart = database2)
#gene2 = getBM(c("ensembl_gene_id", "external_gene_name","go_id","name_1006"), mart = genome2)
#
## Prepare data for Gene Set Analysis
total.genes = colnames(datExpr14_cl) # total genes in your dataset
## Analysis bosTau annotation: GO
database2 <- biomaRt::useMart(biomart="ensembl",
                              dataset="btaurus_gene_ensembl",
                              host="http://www.ensembl.org")
gene2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","go_id","name_1006"),
               mart = database2)
#
dim(gene2); length(unique(gene2$ensembl_gene_id)); length(unique(gene2$go_id))
goName = unique(gene2[,c(3,4)]); goName = goName[order(goName$go_id),]; goName = goName[-1,]
GO = goName$go_id
Name = goName$name_1006
genesGO = unique(subset(gene2,go_id != "")$ensembl_gene_id)
### select non-preserved modules
nonpres_index_b = (which(Zsummary14_b < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
GO_results_b = list()
#
#pdf(paste("GO Enrichment in modules bicor 84bicor_c.pdf"))
for (i in c(1:(length(nonpres_modulenames_b)))){
  module_name = nonpres_modulenames_b[i]
  nopresID_GO = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == module_name)])
  sig.genes = nopresID_GO # total genes in the non-preserved module
  N = length(total.genes[total.genes %in% genesGO])
  S = length(sig.genes[sig.genes %in% genesGO])
  out = data.frame(GO=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric())
  pdf(paste("GO Enrichment in modules bicor",module_name,".pdf"))
  for(j in 1:length(GO)){
    gENEs = subset(gene2, go_id == GO[i])$ensembl_gene_id 
    m = length(total.genes[total.genes %in% gENEs])
    s = length(sig.genes[sig.genes %in% gENEs])
    M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
    Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
    tmp = data.frame(GO = GO[j], Name = Name[j], totalG = m, sigG = s, Pvalue = Pval)
    out = rbind(out,tmp)}
  # select those has 4 more gene in common and pvalue smaller thn 0.05
  ot = subset(out,totalG > 4 & Pvalue < 0.05)
  final = ot[order(ot$Pvalue),];colnames(final) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue")
  GO_results_b[[i]] = final
  print(final %>%
          top_n(dim(final)[1], wt= -pvalue)%>%
          mutate(hitsPerc = Significant_Genes*100/Total_Genes) %>% ## signi genes, v1 = all genes in the go.
          ggplot(aes(x = hitsPerc,
                     y = GO_Name,
                     colour = pvalue,
                     size = Significant_Genes)) +
          #xlim(0,)+
          geom_point() +
          theme_gray()+
          labs(title= paste("GO Enrichment in module",module_name), x="Hits (%)", y="GO term", colour="p value", size="Count")+
          theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
          theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  dev.off()
}
save(GO_results_b, file = "GO_results_b.RData")
print("Step11 - GO finished and data saved")
