# Founction preparation
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA')
source("Function_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net')
load("data_expr_all_with0prepare_no_corrections_top50.RData")
load("permutedStats-actualModules.RData")
load("modulePreservation_methionine.RData")
load("SoftThres_control.RData")
load("modulePreservation_methionine.RData")
load("module_colorsNlabels_control.RData")


load("Enrich_Ensentials.RData")
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b
load("Ensembl2Entrez_Convert.RData")
# Run loops
#===========================================================================================
#                             14. Interpro enrichment                                    ##
#===========================================================================================
Interpro_Enrich_Results_thres005_1102 = 
  InterPro_Enrich(total_genes_all = Total_list_out_ens,
                  sig_genes_all = Sig_list_out_ens,
                  TestingSubsetNames = TestingSubsetNames,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_0113")

############################################################
### =======             Interpro               ========== ##
############################################################
load("Interpro_Enrichment_0113.RData")

# get loop index
all_module = character()
for (i in seq_along(names(Interpro_results_b))){
  all_module[i] = unlist(strsplit(names(Interpro_results_b)[i]," "))[1]
}

#names(Interpro_results_b[[1]])
# loop
all_interpro_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(Interpro_results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_interpro_results[[i]] = tmp_results
  names(all_interpro_results)[i] = all_module[i]
}

require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_interpro_results,file = "Interpro_Results_all_0113.xlsx")

