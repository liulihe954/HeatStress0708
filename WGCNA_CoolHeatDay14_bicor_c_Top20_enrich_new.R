#================================================================================================
###                               0. pkg & functions prep                                  ######
#================================================================================================
source("Functions_Source.R") 
#setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20/Day14_test/")
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
# Following are for testing 'wgcna' --- ignore unless you find them relevent
#ENS_ID_all <- colnames(datExpr14_cl)
nonpres_index_b = (which(Zsummary14_b < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]

#KEGG_results_b = list()
#Kegg_Enrichment_Results = Kegg_Enrich_Plot(ENS_ID_all,
#                                           KEGGthres = 0.05,
#                                           TestingGroupAssignment = moduleColors14_b_cl, 
#                                           TestingSubsetNames = nonpres_modulenames_b,
#                                           keyword = "KEGG_Enrichment_Day14_bicor_c_new")
#===========================================================================================
#                             11. Gene Ontology enrichment                                ##
#===========================================================================================
total.genes = colnames(datExpr14_cl)# total genes in your dataset
GO_results_b = list()
#GO_Enrichment_Results = Go_Enrich_Plot(total.genes,
#                                       GOthres = 0.05,
#                                       TestingGroupAssignment = moduleColors14_b_cl,
#                                       TestingSubsetNames = nonpres_modulenames_b,
#                                       keyword = "GO_Enrichment_Day14_bicor_c_new_z005_0920")

#===========================================================================================
#                               12. Interpro enrichment                                   ##
#===========================================================================================
#nonpres_modulenames_b
ch.total.genes = list();np.genes = list()
for (i in seq_along(nonpres_modulenames_b)){
  np.genes[[i]] = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == nonpres_modulenames_b[i])])
  ch.total.genes[[i]] = colnames(datExpr14_cl)
  names(np.genes)[i] = nonpres_modulenames_b[i]
  names(ch.total.genes)[i] = nonpres_modulenames_b[i]
}

Interpro_Enrichment_Results= InterPro_Enrich(
  total_genes_all =ch.total.genes,
  sig_genes_all = np.genes ,
  TestingSubsetNames = nonpres_modulenames_b,
  IPthres = 0.05,
  biomart="ensembl",
  dataset="btaurus_gene_ensembl",
  Identifier = "ensembl_gene_id",
  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
  keyword = "Interpro_Enrichment_Day14_bicor_c_005_0927")
 
#########################################################################################################################



