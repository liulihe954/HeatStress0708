# getwd()
#================================================================================================
###                                       0. pkg prep                                      ######
#================================================================================================
require(WGCNA);require(sva)
require(ppcor);require(dplyr)
require(edgeR);require(clusterProfiler)
require(ggplot2);require(magrittr)
require(biomaRt);require(gage);require(doParallel)
require(limma);require(recount);require(pamr)
#library(igraph);;library(ggplot2);library(gdata)
#library(ggpubr);require(cowplot);library(extrafont)
#library(plotly) ;library(geomnet);library(readxl);
#library(car) ## qqplot but didnt use
#library(qqplotr)## qqplot: used 
#================================================================================================
###                                       1. dataprep                                      ######
#================================================================================================
#================================================================================================
###                                       1. dataprep                                      ######
#================================================================================================
load("networkData42prepare with corrections_top20.RData")
load('AllColumnsLables.RData')
networkData14_final <- networkData_final
datExpr14_cl = t(networkData14_final[,c(1:(.5*ncol(networkData14_final)))])
datExpr14_ht = t(networkData14_final[,c(((.5*ncol(networkData14_final)+1)):(ncol(networkData14_final)))])
datExpr14_cl = data.frame(datExpr14_cl);datExpr14_ht = data.frame(datExpr14_ht)
dim(datExpr14_cl);dim(datExpr14_ht)
print("datapre done!")
#================================================================================================
###                                  2. weighted in day 14                                ######    
#================================================================================================
## pick soft thresholds

print("Step2 - soft thre plotted and Rdata saved")
#==============================================================================================
#                                 3. soft-threshold and dissimilarity                    ######
#==============================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# Peaerson Cor

# save the matrix
print("Step3 - adj matrix created and rdata saved")
load("AllMatrixDay42_bicor_c.RData")
#==========================================================================================
#                                    4.  plot trees                                  ######
#==========================================================================================
# cor
print("Step4 - dissmi plottd and rdata saved")
#=========================================================================================
#                                5.cutting and merging                              ######
#=========================================================================================
# set the minimum module size relatively high:
print("Step5 - mergeing finished")
load("CoolHeatDay42_bicor_c_top20.RData")
#=================================================================================================
#                              6. plotting heatmap                                            ###
#=================================================================================================
print("Step6 - heapmap created")
#===============================================================================================
#                           7. plot cross-condition dendrogram                               ###
#===============================================================================================
print("Step7 - cross condition dendrogram created")
#=================================================================================================
#         8. Qualitatively and quantitatively measure network preservation at the module level  ##
#=================================================================================================
######### To quantify this result module preservation statistics ######################
# data pre and check data structure
#save(mp14_b, file = "CoolHeatDay42_modulePreservation_bicor_c_top20.RData")
load("CoolHeatDay42_modulePreservation_bicor_c_top20.RData")
print("Step8 - mp finished and data saved")
################ output - shortest - only p summmary  ######################
################ output - shortest - only p summmary  ######################
# specify the reference and the test networks
#===========================================================================================
#                                9. mp visualization                                      ##
#===========================================================================================
####### ###### ######   present Z summary ###### ###### ###### ###### ###### 
# Let us now visualize the data.
modColors14_b = rownames(Obs.PreservationStats14_b)
moduleSize14_b = Obs.PreservationStats14_b$moduleSize
# we will omit the grey module (background genes) and the gold module (random sample of genes)
selectModules = !(modColors14_b %in% c("grey", "gold"))
# Text labels for points
point.label14_b = modColors14_b[selectModules]
#Composite preservation statistics
medianRank14_b=Obs.PreservationStats14_b$medianRank.pres
Zsummary14_b=Z.PreservationStats14_b$Zsummary.pres
#
print("Step9 - preservation statistics calculated vis and saved")
############################################################################################################################
####### ###### ######   present individual Z  ###### ###### ###### ###### ###### 
# Re-initialize module color labels and sizes
print("Step9 - all_module_preservation_statistics finished and data saved")

#===========================================================================================
#                                10. KEGG enrichment                                      ##
#===========================================================================================
##prepare pathway - - - bos taurus
sdb = kegg.gsets(species = "bta", id.type = "kegg", check.new=FALSE)
kegg.gs = sdb$kg.sets[sdb$sigmet.id]
#length(sdb$kg.sets);str(kegg.gs)
mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                         dataset="btaurus_gene_ensembl",
                         host="http://www.ensembl.org")
# select and plot
nonpres_index_b = (which(Zsummary14_b < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
KEGG_results_b = list()
pdf("KEGG_Enrichment_in_modules_bicor_c_top20_day42.pdf")
for (i in c(1:length(nonpres_modulenames_b))){
  #i = 3
  module_name = nonpres_modulenames_b[i]
  nopresID = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == module_name)])
  annot <- getBM(attributes = c("entrezgene_id"),
                 filters="ensembl_gene_id",
                 values = nopresID,
                 mart = mart )
  enterID = annot$entrezgene[!is.na(annot$entrezgene)]
  enrich <- enrichKEGG(gene = enterID, organism = 'bta', qvalueCutoff = 0.05, pvalueCutoff = 0.05)
  KEGG_results_b[[i]] = enrich
  # massage on the dataset (append two columns)
  overlap1 = sub('/.*', '',enrich$GeneRatio)
  total1 = sub('/.*', '',enrich$BgRatio)
  enrich_add = cbind(data.frame(enrich),total1 = as.numeric(total1),overlap1 = as.numeric(overlap1))
  # plot
  #pdf(paste("KEGG Enrichment in module",module_name,".pdf"))
  #if (dim(enrich_add)[1] == 0){ paste(module_name,"not found")} else {
  #pdf(paste("KEGG Enrichment in module", module_name,".pdf"))
  print(enrich_add %>%
          top_n(dim(enrich_add)[1], wt= -pvalue)%>%
          mutate(hitsPerc=overlap1*100/total1)%>% ## signi genes, v1 = all genes in the go.
          ggplot(aes(x=hitsPerc,
                     y=Description,
                     colour=pvalue,
                     size=overlap1)) +
          xlim(0,max(enrich_add$overlap1)+5)+
          geom_point() +
          theme_gray()+
          labs(title= paste("KEGG Enrichment in module",module_name), x="Hits (%)", y="Kegg term", colour="p value", size="Count")+
          theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
          theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
}
dev.off()
save(KEGG_results_b, file = "KEGG_results_bicor_c_top20_day42.RData")
print("Step10 - KEGG finished and data saved")

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
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
GO_results_b = list()
#
pdf("GO_Enrichment_in_modules_bicor_c_top20_day42.pdf")
for (i in c(1:(length(nonpres_modulenames_b)))){
  module_name = nonpres_modulenames_b[i]
  nopresID_GO = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == module_name)])
  sig.genes = nopresID_GO # total genes in the non-preserved module
  N = length(total.genes[total.genes %in% genesGO])
  S = length(sig.genes[sig.genes %in% genesGO])
  out = data.frame(GO=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric())
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
}
dev.off()
save(GO_results_b, file = "GO_results_bicor_c_top20_day42.RData")
print("Step11 - GO finished and data saved")
