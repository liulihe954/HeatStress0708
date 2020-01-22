# Founction preparation
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt')
source("Functions_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/HeatStress0708/Network_No_Crt/Day84_no_correct')
load("CoolHeatday14 bicor.RData")
load("CoolHeatDay14_modulePreservation bicor.RData")
load("networkData84prepare_no_corrections_top50.RData")

# preservation stats pre
ref=1; test = 2
Z.PreservationStats=mp14_b$preservation$Z[[ref]][[test]]
Zsummary=Z.PreservationStats$Zsummary.pres
# label pre
nonpres_index_b = (which(Zsummary < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
# convert ensembl to entrez(ncbi)
TestingModAssign = moduleColors14_b_cl
table(TestingModAssign)
bg_gene = rownames(networkData_50var_nocrt)
TestingSubsetNames = nonpres_modulenames_b 

Convert = ConvertNformat(bg_gene,
                         TestingSubsetNames,
                         TestingModAssign,
                         keyword = "Ensembl2Entrez_Convert")
load("Ensembl2Entrez_Convert.RData")

save(Sig_list_out_entrez,
     Total_list_out_entrez,
     nonpres_modulenames_b,
     Sig_list_out_ens,
     file = "Enrich_Ensentials.RData")
