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
#                                           keyword = "KEGG_Enrichment_Day42_bicor_c_new")

#===========================================================================================
#                             11. Gene Ontology enrichment                                ##
#===========================================================================================
total.genes = colnames(datExpr14_cl)# total genes in your dataset
#GO_results_b = list()
#GO_Enrichment_Results = Go_Enrich_Plot(total.genes,
#                                       GOthres = 0.05,
#                                       TestingGroupAssignment = moduleColors14_b_cl,
#                                       TestingSubsetNames = nonpres_modulenames_b,
#                                       keyword = "GO_Enrichment_Day42_bicor_c_new_z005_0920")

#===========================================================================================
#                              12. Interpro enrichment                                    ##
#===========================================================================================
ch.total.genes = list();np.genes = list()
for (i in seq_along(nonpres_modulenames_b)){
  np.genes[[i]] = as.vector(colnames(datExpr14_cl)[which(moduleColors14_b_cl == nonpres_modulenames_b[i])])
  ch.total.genes[[i]] = colnames(datExpr14_cl)
  names(np.genes)[i] = nonpres_modulenames_b[i]
  names(ch.total.genes)[i] = nonpres_modulenames_b[i]
}

#Interpro_Enrichment_Results= InterPro_Enrich(
#  total_genes_all =ch.total.genes,
#  sig_genes_all = np.genes ,
#  TestingSubsetNames = nonpres_modulenames_b,
#  IPthres = 0.05,
#  biomart="ensembl",
#  dataset="btaurus_gene_ensembl",
#  Identifier = "ensembl_gene_id",
#  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
#  keyword = "Interpro_Enrichment_Day42_bicor_c_005_0927")
 

#===========================================================================================
#                                13. Transfer Identifier                                  ##
#===========================================================================================
total_genes_all = ch.total.genes
sig_genes_all = np.genes
TestingSubsetNames = nonpres_modulenames_b

# Here we have the data compilation
# List with lenght 5 (one for each, easier for looping)
#    sig_genes_all
#    total_genes_all
#
# Now we convert
# two steps - 
# 1. using alias2Symbol(): get the symbol
# 2, convert using biomart
# Potentially loose some, but that's life
## transform ensemble ID or External Gene name to EntrezID
ensembl_try=useMart("ENSEMBL_MART_ENSEMBL",dataset="btaurus_gene_ensembl")
attributes_try = c("ensembl_gene_id","entrezgene_accession","entrezgene_id") #"entrezgene_accession"
gene_try = getBM(attributes=attributes_try,mart = ensembl_try)
match_source = dplyr::select(gene_try,ensembl_gene_id,entrezgene_id)
#
for (i in seq_along(ch.total.genes)){
  sig_genes_all[[i]] = as.character(na.omit(match_source$entrezgene_id[match_source$ensembl_gene_id %in% np.genes[[i]]]))
  names(sig_genes_all)[i] = names(np.genes)[i]
  #message("we lost",length(np.genes[[i]])-length(sig_genes_all[[i]]))
  total_genes_all[[i]] = as.character(na.omit(match_source$entrezgene_id[match_source$ensembl_gene_id %in% ch.total.genes[[i]]]))
  #message("we lost",length(ch.total.genes[[i]])-length(total_genes_all[[i]]))
  names(total_genes_all)[i] = names(ch.total.genes)[i]
}

#str(total_genes_all)
#str(sig_genes_all)


# print out
#require(openxlsx)
#write.xlsx(Sig_list_out,file = "test_convert_sig.xlsx")
#write.xlsx(Total_list_out,file = "test_convert_total.xlsx")

# Keep only the entrez ID: then we have one vector for each element of the list (some format as always)
Sig_list_out_entrez = sig_genes_all
Total_list_out_entrez = total_genes_all
#for (i in c(1:5)){
#  Sig_list_out_entrez[[i]] = data.frame(Sig_list_out[[i]])$entrezgene_id
#  names(Sig_list_out_entrez)[i] = names(Sig_list_out)[i]
#  Total_list_out_entrez[[i]] = data.frame(Total_list_out[[i]])$entrezgene_id
#  names(Total_list_out_entrez)[i] = names(Total_list_out)[i]
#}
#===========================================================================================
#                                14. Reactome enrichment                                  ##
#===========================================================================================
nonpres_index_b = (which(Zsummary14_b < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]


# Read in database
# lowest_path
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = dplyr::filter(NCBI2Reactome_all_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>%
  dplyr::rename(EntrezID=V1, ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V3,V4,V6,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2, 
                Reaction_Description = V3,
                ProteinID = V4,
                Protein_Description = V6,
                Source = V7, Species = V8)

# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
# data massage done 

###
### just for testing
#Total_list_out_entrez_test = Total_list_out_entrez[1:2]
#Sig_list_out_entrez_test = Sig_list_out_entrez[1:2]
#TestingSubsetNames
#InputSource = NCBI2Reactome_all_react_bt
## testing ends

# here we have out list of genes
Sig_list_out_entrez = sig_genes_all
Total_list_out_entrez = total_genes_all

## all react
Reactome_Enrich_all_react_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = nonpres_modulenames_b,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_1001_day42")
## lowest path
Reactome_Enrich_lowest_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = nonpres_modulenames_b,
                                                   InputSource=  NCBI2Reactome_alowest_path_bt,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_1001_day42")
## all path
Reactome_Enrich_all_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = nonpres_modulenames_b,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_1001_day42")

#===========================================================================================
#                                 15. Mesh enrichment                                   ##
#===========================================================================================
#####################
##  run analysis   ##
####################
# just in case that does not work
keyword = "MESH_Enrichment_1001_day42"
DB = paste(keyword,".RData",sep = "")
load(DB)
###
MESH_Enrich_Result1001 = MESH_Enrich(total_genes_all= Total_list_out_entrez,
                                     sig_genes_all = Sig_list_out_entrez,
                                     TestingSubsetNames = nonpres_modulenames_b,
                                     Meshthres = 0.05,
                                     dataset="MeSH.Bta.eg.db",
                                     keyword = "MESH_Enrichment_1001_day42")
