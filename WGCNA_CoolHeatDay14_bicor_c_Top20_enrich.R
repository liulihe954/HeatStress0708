# getwd()
#================================================================================================
###                                       0. pkg prep                                      ######
#================================================================================================
require(sva);require(WGCNA)
require(ppcor);require(dplyr)
require(edgeR);require(clusterProfiler)
require(ggplot2);require(magrittr)
require(biomaRt);require(gage);require(doParallel)
require(limma);require(recount);require(pamr)
library(stringr)
#================================================================================================
###                                       1. dataprep                                      ######
#================================================================================================
load("networkData14prepare with corrections_top20.RData")
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
# Choose a set of soft-thresholding powers
# Mean connectivity as a function of the soft-thresholding power
#save(sft_b_cl,softPower_b,MeanK_b,file = "SoftThres_bicor_c_day14.RData")
print("Step2 - soft thre plotted and Rdata saved")
#==============================================================================================
#                                 3. soft-threshold and dissimilarity                    ######
#==============================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# save the matrix
#save(adjacency14_b_cl,dissTOM14_b_cl,geneTree14_b_cl,adjacency14_b_ht,dissTOM14_b_ht,geneTree14_b_ht,file = "AllMatrixDay14_bicor_c.RData")
#load("AllMatrixDay14_bicor_c.RData")
print("Step3 - adj matrix created and rdata saved")
#==========================================================================================
#                                    4.  plot trees                                  ######
#==========================================================================================
# cor
print("Step4 - dissmi plottd and rdata saved")
#=========================================================================================
#                                5.cutting and merging                              ######
#=========================================================================================
# set the minimum module size relatively high:
print("Step5 - cutting finished")
### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
# Save module colors and labels for use in subsequent parts
#save(MEs14_b_cl, moduleLabels14_b_cl, moduleColors14_b_cl, geneTree14_b_cl, file = "CoolHeatDay_bicor_c_day14.RData")
load("CoolHeatDay_bicor_c_day14.RData")

print("Step5 - mergeing finished")
#=================================================================================================
#                              6. plotting heatmap                                            ###
#=================================================================================================
#===============================================================================================
#                           7. plot cross-condition dendrogram                               ###
#===============================================================================================
print("Step7 - cross condition dendrogram created")
#=================================================================================================
#         8. Qualitatively and quantitatively measure network preservation at the module level  ##
#=================================================================================================
######### To quantify this result module preservation statistics ######################
# data pre and check data structure
#save(mp14_b, file = "CoolHeatDay14_modulePreservation_bicor_c.RData")
load("CoolHeatDay14_modulePreservation_bicor_c.RData")
#load("CoolHeatDay14_modulePreservation_bicor_c.RData")
# load("CoolHeatDay14_modulePreservation.RData")
print("Step8 - mp finished and data saved")
################ output - shortest - only p summmary  ######################
#write.csv(Results_b_mp14_1,"module_size_and_preservation_statistics_bicor_day14.csv")
################ output - shortest - only p summmary  ######################
# specify the reference and the test networks
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
# Look at the observed preservation statistics
# Obs.PreservationStats
# Z statistics from the permutation test analysis Z.PreservationStats
#write.csv(Obs.PreservationStats14_b,"Obs.PreservationStats_bicor_c_day14.csv")
#write.csv(Z.PreservationStats14_b,"Z.PreservationStats_bicor_c_day14.csv")
print("Step8 - preservation statistics calculated and saved")
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
ref = 1;test = 2
# Module labels and module sizes are also contained in the results
modColors14_b = rownames(statsZ14_b)
moduleSizes14_b = mp14_b$quality$Z[[ref]][[test]][,1];
# Exclude improper modules / leave grey and gold modules out
plotMods = !(modColors14_b %in% c("grey", "gold"));
# Create numeric labels for each module
labs14_b= match(modColors14_b[plotMods], standardColors());

# Compare preservation to quality:
# Text labels for points
#text = modColors14_b[plotMods];
# Auxiliary convenience variable
#plotData_b = cbind(mp14_b$preservation$observed[[ref]][[test]][,2], mp14_b$preservation$Z[[ref]][[test]][,2])
# Main titles for the plot
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
# Plot each Z statistic in a separate plot.
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

pdf("KEGG_Enrichment_in_modules_bicor_c_day14.pdf")
for (i in c(1:length(nonpres_modulenames_b))){
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
save(KEGG_results_b, file = "KEGG_results_bicor_day14.RData")
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

pdf("GO_Enrichment_in_modules_bicor_c_day14.pdf")
for (i in c(1:(length(nonpres_modulenames_b)))){
  module_name = nonpres_modulenames_b[i]
  nopresID_GO = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == module_name)])
  sig.genes = nopresID_GO # total genes in the non-preserved module
  N = length(total.genes[total.genes %in% genesGO])
  S = length(sig.genes[sig.genes %in% genesGO]) #
  out = data.frame(GO=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric())
  for(j in 1:length(GO)){
    gENEs = subset(gene2, go_id == GO[i])$ensembl_gene_id # all gene in target GO
    m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
    s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target GO also in the non-preserved module
    M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
    Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
    tmp = data.frame(GO = GO[j], Name = Name[j], totalG = m, sigG = s, Pvalue = Pval)
    out = rbind(out,tmp)}
  # select those has 4 more gene in common and pvalue smaller than 0.05
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
save(GO_results_b, file = "GO_results_bicro_c_day14.RData")
print("Step11 - GO finished and data saved")
