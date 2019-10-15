# getwd()
#================================================================================================
###                                       0. pkg prep                                      ######
#================================================================================================
source("Functions_Source.R")
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
# Choose a set of soft-thresholding powers
print("Step2 - soft thre plotted and Rdata saved")
#==============================================================================================
#                                 3. soft-threshold and dissimilarity                    ######
#==============================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# save the matrix
#save(adjacency14_b_cl,dissTOM14_b_cl,geneTree14_b_cl,adjacency14_b_ht,dissTOM14_b_ht,geneTree14_b_ht,file = "AllMatrixDay42_bicor_c.RData")
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
# save(MEs14_b_cl, moduleLabels14_b_cl, moduleColors14_b_cl, geneTree14_b_cl, file = "CoolHeatDay42_bicor_c_top20.RData")
load("CoolHeatDay42_bicor_c_top20.RData")
print("Step5 - mergeing finished")
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
#stats14_b = mp14_b$preservation$Z$ref.cl$inColumnsAlsoPresentIn.ht
#Results_b_mp14_1 = stats14_b[order(-stats14_b[,2]),c(1:2)]
#save(mp14_b, file = "CoolHeatDay42_modulePreservation_bicor_c_top20.RData")
load("CoolHeatDay42_modulePreservation_bicor_c_top20.RData")
# load("CoolHeatDay14_modulePreservation.RData")
print("Step8 - mp finished and data saved")
################ output - shortest - only p summmary  ######################
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
print(cbind(statsObs14_b[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ14_b[, c("Zsummary.pres", "Zsummary.qual")], 2)))

# Text labels for points
text = modColors14_b[plotMods];
# Auxiliary convenience variable
plotData_b = cbind(mp14_b$preservation$observed[[ref]][[test]][,2], mp14_b$preservation$Z[[ref]][[test]][,2])
# Main titles for the plot
print("Step9 - all_module_preservation_statistics finished and data saved")

#===========================================================================================
#                                10. Data structure pre                                   ##
#                                  format and concersion to entrez
#===========================================================================================
nonpres_index_b = (which(Zsummary14_b < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
#
TestingModAssign42 = moduleColors14_b_cl
bg_gene42 = colnames(datExpr14_cl)
TestingSubsetNames42 = nonpres_modulenames_b 
Convert = ConvertNformat(bg_gene42,
                         TestingSubsetNames42,
                         TestingModAssign42,
                         keyword = "Ensembl2Entrez_Convert_day42")
load("Ensembl2Entrez_Convert_day42.RData")

#str(Sig_list_out_entrez)
#str(Total_list_out_entrez)
#===========================================================================================
#                                10. KEGG enrichment                                      ##
#===========================================================================================
#nonpres_index_b = (which(Zsummary14_b < 2))
#nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
#nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
str(Sig_list_out_entrez)
str(Total_list_out_entrez)
str(Sig_list_out)
Kegg_Enrichment_pval005_1014 = Kegg_Enrich_Plot(sig_genes_all = Sig_list_out_entrez,
                                                total_genes_all = Total_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames42,
                                                KEGGthres = 0.05, 
                                                species = "bta", 
                                                id.type = "kegg",
                                                Sig_list_out =Sig_list_out,
                                                keyword = "Kegg_Enrichment_pval005_1014_Day42")

#==============================================================================================
#                                      11. Mesh enrichment                                   ##
#==============================================================================================
#TestingSubsetNames
MESH_Enrichment_1014 = MESH_Enrich(total_genes_all = Total_list_out_entrez,
                                   sig_genes_all = Sig_list_out_entrez,
                                   TestingSubsetNames = TestingSubsetNames42,
                                   Meshthres = 0.05,
                                   Sig_list_out = Sig_list_out,
                                   MeshCate = c("D","G"),
                                   dataset="MeSH.Bta.eg.db",
                                   keyword = "MESH_Enrichment_1014_Day42")

#===========================================================================================
#                             12. Reactome  enrichment                                    ##
#===========================================================================================
## all react
Reactome_Enrich_all_react_1014 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames42,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Sig_list_out = Sig_list_out,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_1014_Day42")
## lowest path
Reactome_Enrich_lowest_path_1014 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames42,
                                                   InputSource=  NCBI2Reactome_lowest_path_bt,
                                                   Sig_list_out = Sig_list_out,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_1014_Day42")
## all path
Reactome_Enrich_all_path_1014 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames42,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Sig_list_out = Sig_list_out,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_1011_Day14_Day42")
#===========================================================================================
#                             13. Gene Ontology enrichment                                ##
#===========================================================================================
Enrich_Results_thres005_1014 = Go_Enrich_Plot(total_genes_all = Total_list_out_ens,
                                              sig_genes_all = Sig_list_out_ens,
                                              TestingSubsetNames = TestingSubsetNames42,
                                              GOthres = 0.05,
                                              keyword = "GO_Enrichment_pval005_1014_Day14_Day14")

#===========================================================================================
#                             14. Interpro enrichment                                    ##
#===========================================================================================
Interpro_Enrich_Results_thres005_1014 = 
  InterPro_Enrich(total_genes_all = Total_list_out_ens,
                  sig_genes_all = Sig_list_out_ens,
                  TestingSubsetNames = TestingSubsetNames42,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_thres005_1014_Day42")

