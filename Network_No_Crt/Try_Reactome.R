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

#### Read in database
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Reactome_db/")
# lowest_path
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = 
  dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,
                ReactomeID = V2,
                Reactome_Description = V4, 
                Source = V5, 
                Species = V6)
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = 
  dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V4,V6,V2,V3,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V4, 
                Reactome_Description = V6,
                ProteinID = V2,
                Protein_Description = V3,
                Source = V7, Species = V8)
#head(NCBI2Reactome_all_react_bt,50)

# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net")
### Container pre
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b
# Run loops
load("Ensembl2Entrez_Convert.RData")
# Run loops
#===========================================================================================
#                             12. Reactome  enrichment                                    ##
#===========================================================================================
## all react
Reactome_Enrich_all_react_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Sig_list_out = Sig_list_out,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_0113")
## lowest path
Reactome_Enrich_lowest_path_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames,
                                                   InputSource=  NCBI2Reactome_lowest_path_bt,
                                                   Sig_list_out = Sig_list_out,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_0113")
## all path
Reactome_Enrich_all_path_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Sig_list_out = Sig_list_out,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_0113")




############################################################
### =======                   Reactome         ========== ##
############################################################
# get loop index
load('Reactome_Enrich_all_path_0113.RData')
all_module = character()
for (i in seq_along(names(Reactome_results_b))){
  all_module[i] = unlist(strsplit(names(Reactome_results_b)[i]," "))[1]
}

all_data =c("Reactome_Enrich_all_path_0113.RData",
            "Reactome_Enrich_lowest_path_0113.RData",
            "Reactome_Enrichment_all_react_0113.RData")
all_keywords=c("Reactome_all_path_Results_all_0113.xlsx",
               "Reactome_lowest_path_Results_all_0113.xlsx",
               "Reactome_all_react_Results_all_0113.xlsx")



for (m in c(1:3)){
  dataset = all_data[m]
  keyword = all_keywords[m]
  load(dataset)
  # loop
  all_r_a_path_results = list()
  for (i in seq_along(all_module)){
    tmp_name = all_module[i]
    tmp_results = Parse_Results(Reactome_results_b[i], keyword= "-")
    if (!(dim(tmp_results)[1] == 0)){
      tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
    }
    all_r_a_path_results[[i]] = tmp_results
    names(all_r_a_path_results)[i] = all_module[i]
  }
  require(openxlsx)
  setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
  write.xlsx(all_r_a_path_results,file=keyword)
  #setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA")
}

