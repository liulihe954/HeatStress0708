setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt')
source("Functions_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/Day14_no_correct')
load("CoolHeatday14 bicor.RData")
load("CoolHeatDay14_modulePreservation bicor.RData")
load("networkData14prepare_no_corrections_top50.RData")
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
Kegg_Enrichment_pval005_1102 = Kegg_Enrich_Plot(sig_genes_all = Sig_list_out_entrez,
                                                total_genes_all = Total_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                KEGGthres = 0.05, 
                                                species = "bta", 
                                                id.type = "kegg",
                                                Sig_list_out =Sig_list_out,
                                                keyword = "Kegg_Enrichment_0121")
############################################################
### =======                KEGG                ========== ##
############################################################
load("Kegg_Enrichment_0121.RData")
# get loop index
all_module = character()
for (i in seq_along(names(KEGG_results_b))){
  all_module[i] = unlist(strsplit(names(KEGG_results_b)[i]," "))[1]
}
# loop
all_kegg_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(KEGG_results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_kegg_results[[i]] = tmp_results
  names(all_kegg_results)[i] = all_module[i]
}
require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/enrich_results/Day14")
write.xlsx(all_kegg_results,file = "KEGG_Results_all_0121.xlsx")

