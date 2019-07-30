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
#library(igraph);;library(ggplot2);library(gdata)
#library(ggpubr);require(cowplot);library(extrafont)
#library(plotly) ;library(geomnet);library(readxl);
#library(car) ## qqplot but didnt use
#library(qqplotr)## qqplot: used 
#================================================================================================
###                                       1. dataprep                                      ######
#================================================================================================
## read
#setwd("/Users/liulihe95/Desktop/HeatStress0708");getwd()
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
save(networkData14_normalized_normfactors,networkData14_normalized,file = "networkData14 norm and factors.RData")

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


# step 7 check for na (may not necessary)
table(is.na(datExpr14_cl));table(is.na(datExpr14_ht))

#================================================================================================
###                                  2. weighted in day 14                                ######    
#================================================================================================
## pick soft thresholds
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft_b_cl = pickSoftThreshold(datExpr14_cl, powerVector = powers, corFnc = "bicor",verbose = 0)
### 10 works good. 10 - 0.832 and corresponging mean connectivity
softPower_b = min(sft_b_cl$fitIndices[,1][which(sft_b_cl$fitIndices[,2] > 0.8)])
# pre_checked
softPower_b = 22
MeanK_b = sft_b_cl$fitIndices[softPower_b,5]
# Plot the results of threshold picking:
sizeGrWindow(9,5)
cex1 = 0.9
pdf(file = "soft_b threshold picking day14.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
     xlab="Soft Threshold bicor (power) 14_cl",ylab="Scale Free Topology Model Fit bicor, unsigned R^2 14_cl",type="n",
     main = paste("Scale independence bicor 14_cl"));
text(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5],
     xlab="Soft Threshold bicor (power) 14_cl",ylab="Mean Connectivity bicor 14_cl",type="n",
     main = paste("Mean connectivity bicor 14_cl"))
text(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=MeanK_b,col="red")
dev.off()
#save
save(sft_b_cl,softPower_b,MeanK_b,file = "SoftThres bicor.RData")
print("Step2 - soft thre plotted and Rdata saved")
#==============================================================================================
#                                 3. soft-threshold and dissimilarity                    ######
#==============================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# Peaerson Cor
adjacency14_b_cl = adjacency(datExpr14_cl,power=softPower_b,type="unsigned",corFnc = "bicor");
diag(adjacency14_b_cl)=0
dissTOM14_b_cl = 1-TOMsimilarity(adjacency14_b_cl, TOMType="unsigned")
geneTree14_b_cl = hclust(as.dist(dissTOM14_b_cl), method ="average")
#
adjacency14_b_ht = adjacency(datExpr14_ht,power=softPower_b,type="unsigned",corFnc = "bicor");
diag(adjacency14_b_ht)=0
dissTOM14_b_ht = 1-TOMsimilarity(adjacency14_b_ht, TOMType="unsigned")
geneTree14_b_ht = hclust(as.dist(dissTOM14_b_ht), method="average")
# save the matrix
save(adjacency14_b_cl,dissTOM14_b_cl,geneTree14_b_cl,adjacency14_b_ht,dissTOM14_b_ht,geneTree14_b_ht,file = "AllMatrixday14 bicor.RData")
print("Step3 - adj matrix created and rdata saved")
#==========================================================================================
#                                    4.  plot trees                                  ######
#==========================================================================================
# cor
pdf("dendrogram14 bicor.pdf",     height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree14_b_cl,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity bicor (14_cl)",labels=FALSE,hang=0.04);
plot(geneTree14_b_ht,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity bicor (14_ht)",labels=FALSE,hang=0.04);
dev.off()
print("Step4 - dissmi plottd and rdata saved")
#=========================================================================================
#                                5.cutting and merging                              ######
#=========================================================================================
# set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods14_b_cl = cutreeDynamic(dendro = geneTree14_b_cl, distM = dissTOM14_b_cl,
                                   cutHeight=0.995, deepSplit = 1, pamRespectsDendro = FALSE,
                                   minClusterSize = minModuleSize);
table(dynamicMods14_b_cl)
# Convert numeric lables into colors
dynamicColors14_b_cl = labels2colors(dynamicMods14_b_cl)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,16)
pdf("dendrogram_bicor_14cl_nomg.pdf",height=8,width=16)
plotDendroAndColors(geneTree14_b_cl, dynamicColors14_b_cl, "Dynamic Tree Cut bicor 14_cl nomg",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors bicor 14_cl nomg")
dev.off()
print("Step5 - cutting finished")
### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList14_b_cl = moduleEigengenes(datExpr14_cl,colors = dynamicColors14_b_cl)
MEs14_b_cl = MEList14_b_cl$eigengenes
#greyMEName = paste(moduleColor.getMEprefix(), "grey", sep = "") 
#if (greyMEName %in% colnames(MEList14_b_cl$eigengenes))  { print("grey found")
#  MEs14_b_cl = removeGreyME(MEList14_b_cl$eigengenes)}
# Calculate dissimilarity of module eigengenes
MEDiss14_b_cl = 1 - cor(MEs14_b_cl)
# Cluster module eigengenes
METree14_b_cl = hclust(as.dist(MEDiss14_b_cl), method = "average");
# Plot the result of module eigengenes
sizeGrWindow(8,16)
pdf("Clustering of module eigengenes bicor 14_cl.pdf",height=8,width=16)
plot(METree14_b_cl, main = "Clustering of module eigengenes bicor 14_cl",
     xlab = "", sub = "")
## We choose a height cut of 0.2, corresponding to correlation of 0.80, to merge
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge14_b_cl = mergeCloseModules(datExpr14_cl, dynamicColors14_b_cl, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors14_b_cl = merge14_b_cl$colors
# Eigengenes of the new merged modules:
mergedMEs14_b_cl = merge14_b_cl$newMEs;
# Rename to moduleColors
moduleColors14_b_cl = mergedColors14_b_cl
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors())
moduleLabels14_b_cl = match(moduleColors14_b_cl, colorOrder) -1;
MEs14_b_cl = mergedMEs14_b_cl;
# Save module colors and labels for use in subsequent parts
save(MEs14_b_cl, moduleLabels14_b_cl, moduleColors14_b_cl, geneTree14_b_cl, file = "CoolHeatday14 bicor.RData")
print("Step5 - mergeing finished")
#=================================================================================================
#                              6. plotting heatmap                                            ###
#=================================================================================================
pdf("Network heatmap plot bicor, all genes 14_cl.pdf",height=8,width=16)
for (i in c(6:12)){
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM_b = dissTOM14_b_cl^i
  # Set diagonal to NA for a nicer plot
  diag(plotTOM_b) = NA
  # Call the plot function
  TOMplot(plotTOM_b, geneTree14_b_cl, moduleColors14_b_cl, main = "Network heatmap plot bicor, all genes")
}
dev.off()
print("Step6 - heapmap created")
#===============================================================================================
#                           7. plot cross-condition dendrogram                               ###
#===============================================================================================
pdf("Gene dendrogram cross condition bicor.pdf",height=8,width=16)
plotDendroAndColors(geneTree14_b_cl, moduleLabels14_b_cl, "Modules", dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors bicor (14_cl)")
plotDendroAndColors(geneTree14_b_ht, moduleLabels14_b_cl, "Modules", dendroLabels=F,hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors bicor (14_ht)")
dev.off()
print("Step7 - cross condition dendrogram created")
#=================================================================================================
#         8. Qualitatively and quantitatively measure network preservation at the module level  ##
#=================================================================================================
######### To quantify this result module preservation statistics ######################
# data pre and check data structure
dim(datExpr14_ht);dim(datExpr14_cl)
setLabels = c("cl", "ht")
multiExpr14_b = list(cl=list(data = adjacency14_b_cl),
                     ht=list(data = adjacency14_b_ht))
multiColor14_b = list(cl = moduleColors14_b_cl)
names(multiExpr14_b) = setLabels

# permutation
mp14_b = modulePreservation(multiExpr14_b,multiColor14_b,referenceNetworks=1,verbose=3,
                            corFnc = "bicor",
                            networkType="unsigned", nPermutations=1000,
                            maxGoldModuleSize = 500, maxModuleSize = 500,
                            calculateQvalue = T,
                            calculateCor.kIMall = T,
                            calculateClusterCoeff = T,
                            indent = 3)

stats14_b = mp14_b$preservation$Z$ref.cl$inColumnsAlsoPresentIn.ht
Results_b_mp14_1 = stats14_b[order(-stats14_b[,2]),c(1:2)]
save(mp14_b, file = "CoolHeatDay14_modulePreservation bicor.RData")
# load("CoolHeatDay14_modulePreservation.RData")
print("Step8 - mp finished and data saved")
################ output - shortest - only p summmary  ######################
write.csv(Results_b_mp14_1,"module size and preservation statistics bicor day14.csv")
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
write.csv(Obs.PreservationStats14_b,"Obs.PreservationStats bicor.csv")
write.csv(Z.PreservationStats14_b,"Z.PreservationStats bicor.csv")
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
pdf("medianRank_Zsummary versus module size bicor.pdf",height = 8, width = 16)
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize14_b[selectModules],medianRank14_b[selectModules],col=1,
     bg=modColors14_b[selectModules], pch = 21,main="medianRank Preservation day14 bicor",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize14_b[selectModules],medianRank14_b[selectModules],point.label14_b,cex=1,offs=0.03)
# plot Zsummary versus module size
plot(moduleSize14_b[selectModules],Zsummary14_b[selectModules], col = 1,
     bg=modColors14_b[selectModules],pch = 21,main="Zsummary Preservation day14 bicor",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize14_b[selectModules],Zsummary14_b[selectModules],point.label14_b,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()
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
print(cbind(statsObs14_b[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ14_b[, c("Zsummary.pres", "Zsummary.qual")], 2)))

# Text labels for points
text = modColors14_b[plotMods];
# Auxiliary convenience variable
plotData_b = cbind(mp14_b$preservation$observed[[ref]][[test]][,2], mp14_b$preservation$Z[[ref]][[test]][,2])
# Main titles for the plot
mains = c("Preservation Median rank bicor", "Preservation Zsummary bicor")

# Start the plot: open a suitably sized graphical window and set sectioning and margins.
# Plot each Z statistic in a separate plot.
pdf("all module preservation statistics nested bicor.pdf",height = 8, width = 16)
par(mfrow = c(4,5))
for (s in 1:ncol(statsZ14_b)){
  min = min(statsZ14_b[plotMods, s], na.rm = TRUE)
  max = max(statsZ14_b[plotMods, s], na.rm = TRUE)
  if (min > -max/5) min = -max/5
  plot( moduleSizes14_b[plotMods], statsZ14_b[plotMods,s], col = 1, bg = modColors14_b[plotMods], pch = 21,
        main = colnames(statsZ14_b)[s],
        cex = 1.7,
        ylab = colnames(statsZ14_b)[s], xlab = "Module size", log = "x",
        ylim = c( min - 0.1 * (max-min), max + 0.1 * (max-min) ),
        xlim = c(20, 1000))
  labelPoints(moduleSizes14_b[plotMods], statsZ14_b[plotMods,s],labs14_b,cex = 0.7, offs = 0.04)
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()
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
KEGG_results_b = list()
pdf("KEGG Enrichment in modules bicor.pdf")
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
save(KEGG_results_b, file = "KEGG_results_b.RData")
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
GO_results_b = list()
#
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
  pdf(paste("GO Enrichment in modules bicor",module_name,".pdf"))
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
