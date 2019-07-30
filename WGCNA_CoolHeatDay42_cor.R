#================================================================================================
###                                       0. pkg prep                                      ######
#================================================================================================
require(sva);require(WGCNA)
require(ppcor);require(dplyr)
require(edgeR);require(clusterProfiler)
require(ggplot2);require(magrittr)
require(biomaRt);require(gage);require(doParallel)
require(limma);require(recount);require(pamr)
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
networkData14 = networkData42               ###################################
#### substitute dataset for multiple scripts ###################################
################################################################################

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
save(networkData14_normalized_normfactors,networkData14_normalized,file = "networkData42 norm and factors.RData")

# step 3 - select var top 50 & 90% of the counts above 10
# top 50 variance 
networkData14_normalized$variance = apply(networkData14_normalized, 1, var)
sum(networkData14_normalized$variance >= quantile(networkData14_normalized$variance,c(.50)))
networkData14_50var = networkData14_normalized[networkData14_normalized$variance >= quantile(networkData14_normalized$variance,c(.50)), ] #50% most variable genes
networkData14_50var$variance <- NULL
dim(networkData14_50var)
# We suggest removing features whose counts are consistently low 
#(for example, removing all features that have a count of less than say 10 in more than 90% of the samples) 
lowcount_keepIndex14 = (rowSums(networkData14_50var >= 10)) >= round(ncol(networkData14_50var)*.9)
table(lowcount_keepIndex14) # 1437 out again
networkData14_filnal = networkData14_50var[lowcount_keepIndex14,]
dim(networkData14_filnal)
# step 4 - select into two groups
datExpr14_cl = t(networkData14_filnal[,names(networkData14_filnal) %in% names(networkData)[column_42_cl] ])
datExpr14_ht = t(networkData14_filnal[,names(networkData14_filnal) %in% names(networkData)[column_42_ht] ])
datExpr14_cl = data.frame(datExpr14_cl);datExpr14_ht = data.frame(datExpr14_ht)
dim(datExpr14_cl);dim(datExpr14_ht)

# step 5 check for na (may not necessary)
table(is.na(datExpr14_cl));table(is.na(datExpr14_ht))
print("check na results cor")
### summary
# raw 27570 
# remove top 40 counts - 4 /// around (5) (27566)
# remove all zero - around 48k () (22699)
# and normalization (no removal)
# top 90 variance and 90% of sample above 10   - 11349 (remain 11350) and - 83 (11267) 
# rest: 11267
####

###### slective testing  ########
# random select 1500 for testing
#set.seed(316)
#random_label = sample(c(1:ncol(datExpr14_ht)),4000,replace = F)
#datExpr14_cl = datExpr14_cl[,random_label]
#datExpr14_ht = datExpr14_ht[,random_label]
#dim(datExpr14_ht);dim(datExpr14_cl)

#================================================================================================
###                                  2. weighted in day 14                                ######    
#================================================================================================
## good genes --- ok
gsg_all = goodSamplesGenes(networkData14_filnal, verbose = 3)
gsg_all$allOK
## plot sample distance
sample_Tree_ht = hclust(dist(datExpr14_ht), method = "average")
sample_Tree_cl = hclust(dist(datExpr14_cl), method = "average")
# plot samples
pdf(file = "Plots_sampleClustering_day42.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
par(mfrow = c(1,2))
par(cex = 0.6)
plot(sample_Tree_ht, main = "Sample clustering to detect outliers ht42 group",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(sample_Tree_cl, main = "Sample clustering to detect outliers cl42 group",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
print("Step2 - sample tree plotted")

## pick soft thresholds

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft_cl = pickSoftThreshold(datExpr14_cl, powerVector = powers, verbose = 0)
### 10 works good. 10 - 0.832 and corresponging mean connectivity
softPower = min(sft_cl$fitIndices[,1][which(sft_cl$fitIndices[,2] > 0.8)])
# pre_checked
softPower = 18
MeanK = sft_cl$fitIndices[softPower,5]
# Plot the results of threshold picking:
sizeGrWindow(9,5)
cex1 = 0.9
pdf(file = "soft_threshold_picking_day42.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_cl$fitIndices[,1], -sign(sft_cl$fitIndices[,3])*sft_cl$fitIndices[,2],
     xlab="Soft Threshold (power) 42_cl",ylab="Scale Free Topology Model Fit, unsigned R^2 42_cl",type="n",
     main = paste("Scale independence 42_cl"));
text(sft_cl$fitIndices[,1], -sign(sft_cl$fitIndices[,3])*sft_cl$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_cl$fitIndices[,1], sft_cl$fitIndices[,5],
     xlab="Soft Threshold (power) 42_cl",ylab="Mean Connectivity 42_cl",type="n",
     main = paste("Mean connectivity 42_cl"))
text(sft_cl$fitIndices[,1], sft_cl$fitIndices[,5], labels=powers,cex=cex1,col="red")
abline(h=MeanK,col="red")
dev.off()
#save
save(sft_cl,softPower,MeanK,file = "SoftThres_day42.RData")
print("Step2 - soft thre plotted and Rdata saved")

#==============================================================================================
#                                 3. soft-threshold and dissimilarity                    ######
#==============================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# Peaerson Cor
adjacency14_cl = adjacency(datExpr14_cl,power=softPower,type="unsigned",corFnc = "cor");
diag(adjacency14_cl)=0
dissTOM14_cl = 1-TOMsimilarity(adjacency14_cl, TOMType="unsigned")
geneTree14_cl = hclust(as.dist(dissTOM14_cl), method ="average")
#
adjacency14_ht = adjacency(datExpr14_ht,power=softPower,type="unsigned",corFnc = "cor");
diag(adjacency14_ht)=0
dissTOM14_ht = 1-TOMsimilarity(adjacency14_ht, TOMType="unsigned")
geneTree14_ht = hclust(as.dist(dissTOM14_ht), method="average")
# save the matrix
save(adjacency14_cl,dissTOM14_cl,geneTree14_cl,adjacency14_ht,dissTOM14_ht,geneTree14_ht,file = "AllMatrixday42.RData")
print("Step3 - adj matrix created and rdata saved")

#==========================================================================================
#                                    4.  plot trees                                  ######
#==========================================================================================
# cor
pdf("dendrogram42.pdf",     height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree14_cl,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (42_cl)",labels=FALSE,hang=0.04);
plot(geneTree14_ht,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (42_ht)",labels=FALSE,hang=0.04);
dev.off()
print("Step4 - dissmi plottd and rdata saved")
#=========================================================================================
#                                5.cutting and merging                              ######
#=========================================================================================
# set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods14_cl = cutreeDynamic(dendro = geneTree14_cl, distM = dissTOM14_cl,
                                 cutHeight=0.995, deepSplit = 1, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize);
table(dynamicMods14_cl)
# Convert numeric lables into colors
dynamicColors14_cl = labels2colors(dynamicMods14_cl)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,16)
pdf("dendrogram_42cl_nomg.pdf",height=8,width=16)
plotDendroAndColors(geneTree14_cl, dynamicColors14_cl, "Dynamic Tree Cut 42_cl nomg",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors 42_cl nomg")
dev.off()
print("Step5 - cutting finished")

### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList14_cl = moduleEigengenes(datExpr14_cl,colors = dynamicColors14_cl)
MEs14_cl = MEList14_cl$eigengenes
#greyMEName = paste(moduleColor.getMEprefix(), "grey", sep = "") 
#if (greyMEName %in% colnames(MEList14_cl$eigengenes))  { print("grey found")
#  MEs14_cl = removeGreyME(MEList14_cl$eigengenes)}
# Calculate dissimilarity of module eigengenes
MEDiss14_cl = 1 - cor(MEs14_cl)
# Cluster module eigengenes
METree14_cl = hclust(as.dist(MEDiss14_cl), method = "average");
# Plot the result of module eigengenes
sizeGrWindow(8,16)
pdf("Clustering_of_module_eigengenes_42_cl.pdf",height=8,width=16)
plot(METree14_cl, main = "Clustering of module eigengenes 42_cl",
     xlab = "", sub = "")
## We choose a height cut of 0.2, corresponding to correlation of 0.80, to merge
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()
print("Step5 - Clustering of module eigengenes plotted")
# Call an automatic merging function
merge14_cl = mergeCloseModules(datExpr14_cl, dynamicColors14_cl, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors14_cl = merge14_cl$colors
# Eigengenes of the new merged modules:
mergedMEs14_cl = merge14_cl$newMEs;
# Rename to moduleColors
moduleColors14_cl = mergedColors14_cl
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors())
moduleLabels14_cl = match(moduleColors14_cl,colorOrder) -1;
table(moduleLabels14_cl)
MEs14_cl = mergedMEs14_cl;
# Save module colors and labels for use in subsequent parts
save(MEs14_cl, moduleLabels14_cl, moduleColors14_cl, geneTree14_cl, file = "CoolHeatday42.RData")
print("Step5 - mergeing finished")

#=================================================================================================
#                              6. plotting heatmap                                            ###
#=================================================================================================
pdf("Network_heatmap_plot_all_genes_42_cl.pdf",height=8,width=16)
for (i in c(6:12)){
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM14_cl^i
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  # Call the plot function
  TOMplot(plotTOM, geneTree14_cl, moduleColors14_cl, main = "Network heatmap plot, all genes")
}
dev.off()
print("Step6 - heapmap created")
#===============================================================================================
#                           7. plot cross-condition dendrogram                               ###
#===============================================================================================
pdf("Gene_dendrogram_cross_condition_day42.pdf",height=8,width=16)
plotDendroAndColors(geneTree14_cl, moduleLabels14_cl, "Modules", dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (42_cl)")
plotDendroAndColors(geneTree14_ht, moduleLabels14_cl, "Modules", dendroLabels=F,hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (42_ht)")
dev.off()
print("Step7 - cross condition dendrogram created")
#=================================================================================================
#         8. Qualitatively and quantitatively measure network preservation at the module level  ##
#=================================================================================================
######### To quantify this result module preservation statistics ######################
# data pre and check data structure
dim(datExpr14_ht);dim(datExpr14_cl)
setLabels = c("cl", "ht")
multiExpr14 = list(cl=list(data = adjacency14_cl),
                   ht=list(data = adjacency14_ht))
multiColor14 = list(cl = moduleColors14_cl)
names(multiExpr14) = setLabels



# permutation
mp14 = modulePreservation(multiExpr14,multiColor14,referenceNetworks=1,verbose=3,
                          networkType="unsigned", nPermutations=1000, 
                          maxGoldModuleSize=500, maxModuleSize=500, 
                          interpolationPlotFile = "CoolHeatDay42_modulePreservation.pdf")

stats14 = mp14$preservation$Z$ref.cl$inColumnsAlsoPresentIn.ht
Results_mp14_1 = stats14[order(-stats14[,2]),c(1:2)]
save(mp14, file = "CoolHeatDay42_modulePreservation.RData")
print("Step8 - mp finished and data saved")
# load("CoolHeatDay14_modulePreservation.RData")
################ output - shortest - only p summmary  ######################
write.csv(Results_mp14_1,"module_size_and_preservation_statistics_day42.csv")
################ output - shortest - only p summmary  ######################
# specify the reference and the test networks
ref=1; test = 2
### print results - short version
statsObs14 = cbind(mp14$quality$observed[[ref]][[test]][,-1], mp14$preservation$observed[[ref]][[test]][,-1])
statsZ14 = cbind(mp14$quality$Z[[ref]][[test]][, -1], mp14$preservation$Z[[ref]][[test]][,-1]);
# Compare preservation to quality:
print(cbind(statsObs14[,c("medianRank.pres", "medianRank.qual")],
            signif(statsZ14[,c("Zsummary.pres", "Zsummary.qual")],2)))

### print results - full
Obs.PreservationStats14= mp14$preservation$observed[[ref]][[test]]
Z.PreservationStats14=mp14$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
# Obs.PreservationStats
# Z statistics from the permutation test analysis Z.PreservationStats
write.csv(Obs.PreservationStats14,"Obs.PreservationStats_day42.csv")
write.csv(Z.PreservationStats14,"Z.PreservationStats_day42.csv")
print("Step8 - preservation statistics calculated and saved")
#===========================================================================================
#                                9. mp visualization                                      ##
#===========================================================================================
####### ###### ######   present Z summary ###### ###### ###### ###### ###### 
# Let us now visualize the data.
modColors14 = rownames(Obs.PreservationStats14)
moduleSize14 = Obs.PreservationStats14$moduleSize
# we will omit the grey module (background genes) and the gold module (random sample of genes)
selectModules = !(modColors14 %in% c("grey", "gold"))
# Text labels for points
point.label14 = modColors14[selectModules]
#Composite preservation statistics
medianRank14=Obs.PreservationStats14$medianRank.pres
Zsummary14=Z.PreservationStats14$Zsummary.pres
#
pdf("medianRank_Zsummary_versus_module_size_day42.pdf",height = 8, width = 16)
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize14[selectModules],medianRank14[selectModules],col=1,
     bg=modColors14[selectModules], pch = 21,main="medianRank Preservation day14",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize14[selectModules],medianRank14[selectModules],point.label14,cex=1,offs=0.03)
# plot Zsummary versus module size
plot(moduleSize14[selectModules],Zsummary14[selectModules], col = 1,
     bg=modColors14[selectModules],pch = 21,main="Zsummary Preservation day14",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize14[selectModules],Zsummary14[selectModules],point.label14,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()
print("Step9 - preservation statistics calculated vis and saved")

############################################################################################################################
####### ###### ######   present individual Z  ###### ###### ###### ###### ###### 
# Re-initialize module color labels and sizes
ref = 1
test = 2
# Module labels and module sizes are also contained in the results
modColors14 = rownames(statsZ14)
moduleSizes14 = mp14$quality$Z[[ref]][[test]][,1];
# Exclude improper modules / leave grey and gold modules out
plotMods = !(modColors14 %in% c("grey", "gold"));
# Create numeric labels for each module
labs14= match(modColors14[plotMods], standardColors());
# Compare preservation to quality:
print(cbind(statsObs14[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ14[, c("Zsummary.pres", "Zsummary.qual")], 2)))
# Text labels for points
text = modColors14[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp14$preservation$observed[[ref]][[test]][,2], mp14$preservation$Z[[ref]][[test]][,2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary")
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
# Plot each Z statistic in a separate plot.
pdf("all_module_preservation_statistics_nested_day42.pdf",height = 8, width = 16)
par(mfrow = c(4,5))
for (s in 1:ncol(statsZ14)){
  min = min(statsZ14[plotMods, s], na.rm = TRUE)
  max = max(statsZ14[plotMods, s], na.rm = TRUE)
  if (min > -max/5) min = -max/5
  plot( moduleSizes14[plotMods], statsZ14[plotMods,s], col = 1, bg = modColors14[plotMods], pch = 21,
        main = colnames(statsZ14)[s],
        cex = 1.7,
        ylab = colnames(statsZ14)[s], xlab = "Module size", log = "x",
        ylim = c( min - 0.1 * (max-min), max + 0.1 * (max-min) ),
        xlim = c(20, 1000))
  labelPoints(moduleSizes14[plotMods], statsZ14[plotMods,s],labs14,cex = 0.7, offs = 0.04)
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
nonpres_index = (which(Zsummary14 < 2))
nonpres_modulenames = rownames(Z.PreservationStats14)[nonpres_index]
KEGG_results = list()
pdf("KEGG_Enrichment_in_modules_day42.pdf")
for (i in c(1:length(nonpres_modulenames))){
  #i = 3
  module_name = nonpres_modulenames[i]
  nopresID = as.vector(colnames(datExpr14_cl)[which(moduleColors14_cl == module_name)])
  annot <- getBM(attributes = c("entrezgene_id"),
                 filters="ensembl_gene_id",
                 values = nopresID,
                 mart = mart )
  enterID = annot$entrezgene[!is.na(annot$entrezgene)]
  enrich <- enrichKEGG(gene = enterID, organism = 'bta', qvalueCutoff = 0.05, pvalueCutoff = 0.05)
  KEGG_results[[i]] = enrich
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
save(KEGG_results, file = "KEGG_results_day42.RData")
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
nonpres_index = (which(Zsummary14 < 2))
nonpres_modulenames = rownames(Z.PreservationStats14)[nonpres_index]
GO_results = list()
#

for (i in c(1:(length(nonpres_modulenames)))){
  module_name = nonpres_modulenames[i]
  nopresID_GO = as.vector(colnames(datExpr14_cl)[which(moduleColors14_cl == module_name)])
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
  GO_results[[i]] = final
  pdf(paste("GO_Enrichment_in_modules",module_name,".pdf"))
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
save(GO_results, file = "GO_results_day42.RData")
print("Step11 - GO finished and data saved")


