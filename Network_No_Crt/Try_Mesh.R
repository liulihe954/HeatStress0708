# Founction preparation
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt')
source("Functions_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/Day84_no_correct')
load("CoolHeatday14 bicor.RData")
load("CoolHeatDay14_modulePreservation bicor.RData")
load("networkData84prepare_no_corrections_top50.RData")
load("permutedStats-actualModules.RData")
load("SoftThres bicor.RData")
#load("module_colorsNlabels_control.RData")

load("Enrich_Ensentials.RData")
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b

load("Ensembl2Entrez_Convert.RData")
# Run loops

# MeSH db pre
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Mesh_db")
keyword_outer = "MeshDB2"
DB = paste(keyword_outer,".RData",sep = "")
load(DB)
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/Day14_no_correct')
# Run loops
#==============================================================================================
#                                      11. Mesh enrichment                                   ##
#==============================================================================================
#TestingSubsetNames
MESH_Enrichment_1102 = MESH_Enrich(total_genes_all = Total_list_out_entrez,
                                   sig_genes_all = Sig_list_out_entrez,
                                   TestingSubsetNames = TestingSubsetNames,
                                   Meshthres = 0.05,
                                   Sig_list_out = Sig_list_out,
                                   MeshCate = c("D","G"),
                                   dataset="MeSH.Bta.eg.db",
                                   keyword = "MESH_Enrichment_0121")

############################################################
### =======                Mesh                ========== ##
############################################################
load("MESH_Enrichment_0121.RData")
# get loop index
all_module = character()
for (i in seq_along(names(Mesh_results_b))){
  all_module[i] = unlist(strsplit(names(Mesh_results_b)[i]," "))[1]
}
# loop
all_mesh_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(Mesh_results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_mesh_results[[i]] = tmp_results
  names(all_mesh_results)[i] = all_module[i]
}
require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/enrich_results/Day14")
write.xlsx(all_mesh_results,file = "Mesh_Results_all_0121.xlsx")
