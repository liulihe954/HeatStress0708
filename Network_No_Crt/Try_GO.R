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
#                             13. Gene Ontology enrichment                                ##
#===========================================================================================
Enrich_Results_thres005_1102 = Go_Enrich_Plot(total_genes_all = Total_list_out_ens,
                                              sig_genes_all = Sig_list_out_ens,
                                              TestingSubsetNames = TestingSubsetNames,
                                              GOthres = 0.05,
                                              keyword = "GO_Enrichment_0113")

##################################################
### =======             GO           ========== ##
##################################################
load("GO_Enrichment_0113.RData")
# get loop index

all_module = character()
for (i in seq_along(names(GO_results_b))){
  all_module[i] = unlist(strsplit(names(GO_results_b)[i]," "))[1]
}
# get names
biomart="ensembl";dataset="btaurus_gene_ensembl";attributes = c("go_id","namespace_1003")
database = useMart(biomart);genome = useDataset(dataset, mart = database);gene = getBM(attributes,mart = genome)
namespace_index = dplyr::filter(gene,go_id != "",namespace_1003 != "")
# loop
all_go_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_go_results = Parse_Results(GO_results_b[i], keyword= "-")
  tmp_go_results = dplyr::select(tmp_go_results,-ExternalLoss_total,-InternalLoss_sig) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join( namespace_index, by =c("GOID" ="go_id")) %>% 
    dplyr::rename(go_id = GOID)
  all_go_results[[i]] = tmp_go_results
}
names(all_go_results) = all_module

require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_go_results,file = "GO_Results_all_0113.xlsx")

###
# setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/enrich_results")
# pdf("GO_Enriche_Semantic_Sim.pdf")
# for (m in seq_along(all_go_results)){
#   tmp_out = data.frame(all_go_results[[i]])
#   ReduceDim_GO_Plot(tmp_out,
#                     GOthres = 0.05,
#                     label_sizeCC = 0.4,
#                     label_sizeBP = 0.4,
#                     label_sizeMF = 0.4,
#                     Database = "org.Bt.eg.db",
#                     measure="Jiang",combine=NULL,
#                     Dataset_Name = "Methylation_GO_Enrich")
# }
# dev.off()
# setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA")





