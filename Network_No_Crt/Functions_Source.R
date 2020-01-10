#===========================================================================================
#                                0. Package preparations                                  ##
#===========================================================================================
#devtools::install_github("jakesherman/easypackages")
#my_packages <- c("WGCNA","ppcor","edgeR","clusterProfiler","magrittr","gage","doParallel","recount","pamr","stringr")
#libraries(my_packages)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(limma)
library(sva)
library(recount)
library(biomaRt)
library(readxl)
library(GOSemSim)
library(corrplot)
library(org.Bt.eg.db)
library(meshr)
library(MeSH.db)
library(MeSH.Bta.eg.db)
#===========================================================================================
#                                1. "Massage" - Data processing                           ##
#===========================================================================================
# Classical way, networkData(rawdata,genes in rows and samples in cols), 
#                 n1-n2 (number of ref and treatmnent group)
#                 perct (the percentage to REMOVE based on VARIANCE across two conditions)
# Same rationales but FANCY way, remove confounding artifacts
# (ref - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1700-9 )
DataPre_C = function(networkData, cousin = 0.4, n1, n2, perct,
                     thres_rmzero = 5,count_rmzero,
                     Correct = 'Y'){
  #function prepare
  check_zero = function(networkData,thres_rmzero,count_rmzero){
    cow_count_index = rep("ok",length(rownames(networkData)))
    for (i in seq_along(rownames(networkData))){
      tmp_count = sum(networkData[i,] <= thres_rmzero)
      if (tmp_count >= count_rmzero){cow_count_index[i] = "out"}
    }
    return(cow_count_index)
  }
  #function prepare
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
  q_normalize <- function(dat){
    n = nrow(dat)
    p = ncol(dat)
    rank.dat =  dat # matrix for ranking
    for (i in 1:p){
      rank.dat[,i] = rank(dat[,i])
    }
    U = rank.dat/(n+1)
    qnorm(U)
  }
  Correct_pca = function(rse_raw,method){
    rse_raw <- t(rse_raw)# transpose data so that rows are samples and columns are gene expression measurements
    mod=matrix(1,nrow=dim(rse_raw)[1],ncol=1)
    colnames(mod)="Intercept"
    ## num.sv requires data matrix with features(genes) in the rows and samples in the column
    nsv=num.sv(t(rse_raw), mod, method = method)
    print(paste("Number of PCs estimated to be removed:", nsv))
    ## PC residualization of gene expression data using sva_network. Rows correspond to samples, Columns correspond to features
    exprs_corrected = sva_network(rse_raw, nsv)
    ## Quantile normalize the corrected gene expression measurements such that each expression of every gene follows a gaussian distribution
    exprs_corrected_norm <- q_normalize(exprs_corrected)
    return(list(exprs_corrected_norm = t(data.frame(exprs_corrected_norm))))
  }
  # 1/2/3 normalization then var select for non-correction option
  # 
  # step 1 - filter out top 40% counts
  ## filter out top 40% counts # function established for future use
  networkData_filter = remove_filter(networkData,cousin)$networkData_filter
  # step 2 - normalization (0s out and normalization)
  remove_index = which(rowSums(networkData_filter) == 0)#;length(remove_index14)
  networkData_nm1 = networkData_filter[-remove_index,]#;dim(networkData_nm1)
  networkData_nmList = DGEList(counts = networkData_nm1,group  = c(rep("ref",n1),rep("test",n2)))
  networkData_nm2 = calcNormFactors(networkData_nmList)
  networkData_normalized_normfactors = networkData_nm2$samples
  networkData_normalized = data.frame(networkData_nm2$counts)
  # step 3 - check zeros and var selection for non-correction option
  zero_cm_label = check_zero(networkData_normalized,
                             thres_rmzero = thres_rmzero,
                             count_rmzero = count_rmzero)
  networkData_normalized_nozero = networkData_normalized[zero_cm_label=="ok",]
  #
  networkData_normalized_nozero$variance = apply( networkData_normalized_nozero, 1, var)
  networkData_50var_nocrt =  networkData_normalized_nozero[ networkData_normalized_nozero$variance >= quantile( networkData_normalized_nozero$variance,c(perct)),] #50% most variable genes
  networkData_50var_nocrt$variance <- NULL
  # step 4- log 2 trans
  # log trans
  networkData_log2 = log2(networkData_normalized_nozero+2)
  # step 4 - filter out bottom xx% variation
  # crt
  networkData_log2$variance = apply(networkData_log2,1,var)
  networkData_log2_50var = networkData_log2[networkData_log2$variance >= quantile(networkData_log2$variance,c(perct)),] 
  networkData_log2_50var$variance <- NULL
  # step 5 - pca correction 
  networkData_correction = Correct_pca(networkData_log2_50var,"leek")
  networkData_final = data.frame(networkData_correction$exprs_corrected_norm); 
  names(networkData_final) = names(networkData_log2_50var)
  if (Correct == "Y"){
    save(networkData_final,
         networkData_log2_50var,
         networkData_normalized_normfactors,
         networkData_normalized,
         file = paste(deparse(substitute(networkData)),"prepare_with_corrections","_top",100*(1-perct),".RData",sep = ""))
    return(list(Corrected_log2_PC = networkData_final))}
  else if (Correct == "N"){
    message('pc correction not applied')
    save(networkData_50var_nocrt,
         networkData_normalized_normfactors,
         networkData_normalized,
         file = paste(deparse(substitute(networkData)),"prepare_no_corrections","_top",100*(1-perct),".RData",sep = ""))
    return(list(networkData_50var_no_crt = networkData_50var_nocrt))}
  else {message("please specify pc data correction option - Correct = 'Y' or 'N'")}
}
#########################################################################################################################

#===========================================================================================
#                                2. KEGG enrichment                                      ##
#===========================================================================================
##########################################################################################
### Function necessities: a."ENS_ID_all", which is the collection of all genes in your dataset (format - vector - ensembl iD)
###                      "nonpres_modulenames_b", which is the "sig" module names (format - char - e.g. color names)
###                      "moduleColors14_b_cl", which is the color/module assignments (format - char - same len with ENS_ID_all - this is the module assignment)
###        
###                  Motivation; double loop - 
###                  for (i in seq_along(#modules)){
###                         for (j in seq_along(#all_KEGG))
###                              "fisher.test"$'signiciant_hits'
###                   }
### Function output: a. a plot with all the enriched items, one page for each module.
###                  b. .RData containing a long list, each element show all the enriched items for corresponding module.
###                      same length with that of all non-preserved module
###
##########################################################################################
Parse_Results = function(Results_List,keyword = "Which D.B"){
  all_enrich = data.frame()
  for (i in 1:length(Results_List)){
    len = dim(data.frame(Results_List[i]))[1]
    if (len> 0){
      tmp = data.frame(Results_List[i])
      names(tmp) = names(Results_List[[1]])
      all_enrich = rbind(all_enrich,tmp)
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich)[1]
  total_modules = length(Results_List)
  print(paste("In database: ",keyword,"-",total_hits,"hits found in",total_modules,"tested modules: ",names(Results_List)))
  return(ParseResults = all_enrich)
}

##############################################################################################################
Go_Enrich_Plot = function(total_genes_all,
                          sig_genes_all,
                          TestingSubsetNames,
                          GOthres = 0.05,
                          biomart="ensembl",
                          dataset="btaurus_gene_ensembl",
                          attributes = c("ensembl_gene_id","go_id","name_1006"),
                          keyword = "GO_Enrichment_thres_point1_5sets"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  GO_results_b = list()
  GO_results_b_raw = list()
  DB_List = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr)# load pkg
  ## Analysis bosTau annotation: GO
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  goName = unique(gene[,c(2,3)]); goName = goName[order(goName$go_id),];goName = goName[-1,]
  GO = goName$go_id
  Name = goName$name_1006
  genesGO = unique(subset(gene,go_id != "")$ensembl_gene_id)[-1]#
  for ( p in seq_along(GO)){
    tmp = subset(gene, go_id == GO[p])$ensembl_gene_id
    DB_List[[p]] = tmp #
    names(DB_List)[p]  <- paste(GO[p],"-",Name[p])
  }
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of GO sets to check: ",length(GO)," with total number of names: ",length(Name))
  # plot
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesGO])
    S = length(sig.genes[sig.genes %in% genesGO]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(GO=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(GO)){
      if (j%%100 == 0) {message("tryingd on GO ",j," - ",GO[j]," - ",Name[j])}
      gENEs = DB_List[[j]] # all gene in target GO #### note
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      PastefindG = paste(findG,collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(GO = GO[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    GO_results_b_raw[[i]] = final_raw; names(GO_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched GO raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= GOthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    GO_results_b[[i]] = final;names(GO_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched GO")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_  gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(DB_List,GO_results_b, GO_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant GO terms found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",GOthres)
  message("Nice! - GO enrichment finished and data saved")}
#########################################################################################################################
ReduceDim_GO_Plot = function(Enrich_Out,
                             GOthres = 0.001,
                             label_sizeCC = 0.4,
                             label_sizeBP = 0.4,
                             label_sizeMF = 0.4,
                             Database = "org.Bt.eg.db",
                             measure="Jiang",combine=NULL,
                             Dataset_Name){
  # load libraries + download ref database
  library(GOSemSim);library(corrplot);library(tidyverse)
  do.call(library,list(Database))
  semData_BP <- godata(paste(Database), ont="BP", computeIC=T)
  semData_MF <- godata(paste(Database), ont="MF", computeIC=T)
  semData_CC <- godata(paste(Database), ont="CC", computeIC=T)
  # selection + formating: for each category we have one vector containing all the sig GO terms
  BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>% 
    dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL # name is an attribute and we dont them, so set null
  CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>% 
    dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
  MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>% 
    dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
  ### Now we are trying to get all similarity matrix ready. N x N, symetric, diag = 1
  # For BP
  
  goSimMatrix_BP = GOSemSim::mgoSim(BP_List,
                                    BP_List,
                                    semData=semData_BP,measure=measure,combine = combine)
  suspectID_BP = rownames(goSimMatrix_BP)[is.na(goSimMatrix_BP[,1])]
  if (length(suspectID_BP) != 0){BP_List_new = setdiff(BP_List,suspectID_BP)
  message(length(suspectID_BP)," invalid ID captured in BP: ",suspectID_BP,", thus been removed!")
  } else {BP_List_new = BP_List;message("Nice! All IDs are valid in BP!")}
  goSimMatrix_BP_new = GOSemSim::mgoSim(BP_List_new,
                                        BP_List_new,
                                        semData=semData_BP,measure=measure,combine = combine)
  colnames(goSimMatrix_BP_new) = paste(BP_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)])
  rownames(goSimMatrix_BP_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)],BP_List_new)
  # For CC
  goSimMatrix_CC = GOSemSim::mgoSim(CC_List,
                                    CC_List,
                                    semData=semData_CC,measure=measure,combine = combine)
  suspectID_CC = rownames(goSimMatrix_CC)[is.na(goSimMatrix_CC[,1])]
  if (length(suspectID_CC) != 0){CC_List_new = setdiff(CC_List,suspectID_CC)
  message(length(suspectID_CC)," invalid ID captured in CC: ",suspectID_CC,", thus been removed!")
  } else {CC_List_new = CC_List;message("Nice! All IDs are valid in CC!")}
  goSimMatrix_CC_new = GOSemSim::mgoSim(CC_List_new,
                                        CC_List_new,
                                        semData=semData_CC,measure=measure,combine =combine)
  colnames(goSimMatrix_CC_new) = paste(CC_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)])
  rownames(goSimMatrix_CC_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)],CC_List_new)
  # For MF
  goSimMatrix_MF = GOSemSim::mgoSim(MF_List,
                                    MF_List,
                                    semData=semData_MF,measure=measure,combine = combine)
  suspectID_MF = rownames(goSimMatrix_MF)[is.na(goSimMatrix_MF[,1])]
  if (length(suspectID_MF) != 0){MF_List_new = setdiff(MF_List,suspectID_MF)
  message(length(suspectID_MF)," invalid ID captured in MF: ",suspectID_MF,", thus been removed!")
  } else {MF_List_new = MF_List;message("Nice! All IDs are valid in MF!")}
  goSimMatrix_MF_new = GOSemSim::mgoSim(MF_List_new,
                                        MF_List_new,
                                        semData=semData_MF,measure=measure,combine = combine)
  colnames(goSimMatrix_MF_new) = paste(MF_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)])
  rownames(goSimMatrix_MF_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)],MF_List_new)
  # Now we take the results and plot
  pdf(paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".pdf",sep = ""))
  corrplot(goSimMatrix_CC_new,title = "Semantic_Similarity_Measure_CC",
           tl.col = "black", tl.cex = label_sizeCC, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_BP_new,title = "Semantic_Similarity_Measure_BP",
           tl.col = "black", tl.cex = label_sizeBP, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_MF_new,title = "Semantic_Similarity_Measure_MF",
           tl.col = "black", tl.cex = label_sizeMF, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  dev.off()
  message(dim(goSimMatrix_CC_new)[1],",",
          dim(goSimMatrix_BP_new)[1],",",
          dim(goSimMatrix_MF_new)[1]," GOs ploted in CC, BP and MF, respectively",
          ", cutting at, ",GOthres)
  require(openxlsx)
  CorMatrix <- list("CorMat_BP" = data.frame(goSimMatrix_BP), 
                    "CorMat_CC" = data.frame(goSimMatrix_CC),
                    "CorMat_MF" = data.frame(goSimMatrix_MF))
  write.xlsx(CorMatrix, row.names=TRUE,
             file = paste("Semantic_Similarity_Measure_",
                          Dataset_Name,"_",
                          formatC(GOthres, format = "e", digits = 0),".xlsx",sep = ""))
  save(goSimMatrix_CC,goSimMatrix_BP,goSimMatrix_MF,
       file = paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".RData",sep = ""))
  message("Nice! Excels, Plots exported and RData saved!")
}
#########################################################################################################################
InterPro_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           IPthres = 0.05,
                           biomart="ensembl",
                           dataset="btaurus_gene_ensembl",
                           Identifier = "external_gene_name",
                           attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Interpro_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Interpro_results_b = list()
  Interpro_results_b_raw = list()
  DB_List = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg
  ## GetInterpro : bosTau 
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  ##
  InterproName = unique(gene[,c("interpro","interpro_description")]) %>% arrange(interpro)
  Interpro = na.omit(InterproName$interpro)[-1]
  Name = na.omit(InterproName$interpro_description)[-1]
  #
  if (Identifier == "ensembl_gene_id"){
    genesInterpro = unique(subset(gene,interpro != "")$ensembl_gene_id)
    for ( p in seq_along(Interpro)){
      tmp = subset(gene, interpro == Interpro[p])$ensembl_gene_id
      DB_List[[p]] = tmp #
      names(DB_List)[p]  <- paste(Interpro[p],"-",Name[p])
    }
  } else if (Identifier == "external_gene_name") {
    genesInterpro= unique(subset(gene,interpro != "")$external_gene_name);
    genesInterpro = genesInterpro[-1]
    for ( p in seq_along(Interpro)){
      tmp = subset(gene, interpro == Interpro[p])$external_gene_name
      DB_List[[p]] = tmp #
      names(DB_List)[p]  <- paste(Interpro[p],"-",Name[p])
    }
  } else {message("Sorry, we only have ensembel and names available as identifier, please use one of the followings: 
                  ensembl_gene_id OR external_gene_name.")}
  
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Interpro domains to check: ",length(Interpro)," with total number of names: ",length(Name))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesInterpro])
    S = length(sig.genes[sig.genes %in% genesInterpro]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(Interpro=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(Interpro)){
      if (j%%100 == 0) {message("tryingd on Interpro ",j," - ",Interpro[j]," - ",Name[j])}
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG) # # genes from target interpro also in the non-preserved module
      PastefindG = paste(findG, collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Interpro = Interpro[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Interpro_results_b_raw[[i]] = final_raw; names(Interpro_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Interpro raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= IPthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Interpro_results_b[[i]] = final;names(Interpro_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Interpro")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Interpro = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Interpro_results_b, Interpro_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Interpro domains found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",IPthres)
  message("Nice! - Interpro enrichment finished and data saved")}
#########################################################################################################################
MESH_Enrich = function(total_genes_all,
                       sig_genes_all,
                       TestingSubsetNames,
                       Meshthres = 0.05,
                       Sig_list_out,
                       MeshCate = c("D","G"),
                       #biomart="ensembl",
                       dataset="MeSH.Bta.eg.db",
                       #dataset= "btaurus_gene_ensembl",
                       #Identifier = "external_gene_name",
                       #attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                       keyword = "MESH_Enrichment"){
  #total.genes = Total_list_out_entrez_test
  #sig.genes = Sig_list_out_entrez_test
  #TestingSubsetNames = TestingSubsetNames_test
  total_enrich = 0                        
  raw_pvalue_all = numeric()
  Mesh_results_b = list()
  Mesh_results_b_raw = list()
  DB_List = list()
  library(MeSH.db);library(MeSH.Bta.eg.db);library(tidyverse);library(gage);library(magrittr)
  library(ggplot2);library(biomaRt);library(MeSH.Bta.eg.db)
  ### Three ways to get meshdb
  # 1 download from github: we are gonna use
  #githubURL <- "https://github.com/liulihe954/Repro_Estrous_0918/raw/master/MeshDB.RData"
  #githubURL <- "https://github.com/liulihe954/HeatStress0708/raw/master/MeshDB_new.RData"
  #load(url(githubURL))
  #if (all(MeshCate%in%c("G","D"))){list_Bta = dplyr::filter(list_Bta, MESHCATEGORY %in% MeshCate)}
  #else {
  #  message("Sorry, we only have G and D")
  #  message("Now reload the category you need, it will take a while...")
  #  KEY = keys(MeSH.db, keytype = "MESHID")
  #  List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  #  Match_List = dplyr::select(List, MESHID, MESHTERM)
  #  key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #  list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
  #    dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% MeshCate) %>% 
  #    dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))}
  # 2. match from the very begining (will take an hour or so)
  #KEY = keys(MeSH.db, keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY[1:3], columns = columns(MeSH.db), keytype = "MESHID")
  #Match_List = dplyr::select(List, MESHID, MESHTERM)
  ##head(Match_List) 
  ### Prepare Bta database
  #key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
  #  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% MeshCate) %>% 
  #  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
  
  # 3. alternatively, if you have them in your environment
  keyword_outer = "MeshDB"
  DB = paste(keyword_outer,".RData",sep = "")
  load(DB)
  #Sig_list_out_entrez_test2
  #Total_list_out_entrez_test2
  # Get index
  list_Bta = list_Bta[which(list_Bta$MESHCATEGORY %in% MeshCate),]
  #list_Bta = dplyr::filter(list_Bta,MESHCATEGORY%in%MeshCate)
  genesMesh = unique(list_Bta$GENEID)
  MeshRecords = unique(list_Bta[,c("MESHID","MESHTERM")]) %>% arrange(MESHID)
  MeshID = na.omit(MeshRecords$MESHID)
  MeshTerm = na.omit(MeshRecords$MESHTERM)
  for ( p in seq_along(MeshID)){
    tmp = subset(list_Bta, MESHID == MESHID[p])$GENEID
    DB_List[[p]] = tmp #
    names(DB_List)[p]  <- paste(MeshID[p],"-",MeshTerm[p])
  }
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Mesh to check: ",length(MeshID)," with total number of names: ",length(MeshTerm))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesMesh])
    S = length(sig.genes[sig.genes %in% genesMesh]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(MeshID=character(),
                     MeshTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in c(1:length(MeshID))){
      if (j%%100 == 0) {message("tryingd on MeshID ",j," - ",MeshID[j]," - ",MeshTerm[j])}
      #target = MeshID[j]
      #gENEs = unique(subset(list_Bta, MESHID == target)$GENEID)
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs]) # genes from target  and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      #length(gENEs);Pval
      tmp = data.frame(MeshID= MeshID[j], 
                       MeshTerm = MeshTerm[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Mesh_results_b_raw[[i]] = final_raw; names(Mesh_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Mesh raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Meshthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Mesh_results_b[[i]] = final;names(Mesh_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Mesh")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Mesh = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Mesh_results_b, Mesh_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant MeshIDs found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",Meshthres)
  message("Nice! - Mesh enrichment finished and data saved")}
#########################################################################################################################
Reactome_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           InputSource,
                           Sig_list_out,
                           Reacthres = 0.05,
                           #biomart="ensembl",
                           #dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           #attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Reactome_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Reactome_results_b = list()
  Reactome_results_b_raw = list()
  DB_List = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg
  Reactome_gene =   unique(InputSource[,c("EntrezID")])
  ReactomeRecords = dplyr::select(InputSource,ReactomeID,Reactome_Description) %>% dplyr::arrange(ReactomeID) %>% distinct()
  #ReactomeRecords = unique(InputSource[,c("ReactomeID","Reactome_Description")]) %>% arrange(ReactomeID) #
  ReactomeID = na.omit(ReactomeRecords$ReactomeID)
  ReactomeName = na.omit(ReactomeRecords$Reactome_Description)
  for ( p in seq_along(ReactomeID)){
    IDindex = ReactomeID[p]
    tmp = subset(InputSource, ReactomeID == IDindex)$EntrezID
    DB_List[[p]] = tmp #
    names(DB_List)[p]  <- paste(ReactomeID[p],"-",ReactomeName[p])
  }
  #ReactomeID = ReactomeID[1:300]
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Reactome to check: ",length(ReactomeID)," with total number of names: ",length(ReactomeName))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    #i = 1
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    #length(sig.genes)
    #length(total.genes)
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% Reactome_gene])
    S = length(sig.genes[sig.genes %in% Reactome_gene]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ReactomeID=character(),
                     ReactomeTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(ReactomeID)){
      # j = 101
      if (j%%100 == 0) {message("tryingd on Reactome ",j," - ",ReactomeID[j]," - ",ReactomeName[j])}
      #target = ReactomeID[j]
      #gENEs = unique(subset(InputSource, ReactomeID == target)$EntrezID)
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs]) 
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(ReactomeID = ReactomeID[j], 
                       ReactomeName = ReactomeName[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Reactome_results_b_raw[[i]] = final_raw; names(Reactome_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Reactomeid raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Reacthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Reactome_results_b[[i]] = final;names(Reactome_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched ReactomeID")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Reactome = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Reactome_results_b, Reactome_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Reactome domains found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",Reacthres)
  message("Nice! - Reactome enrichment finished and data saved")}
#########################################################################################################################
Kegg_Enrich_Plot = function(sig_genes_all,
                            total_genes_all,# all genes in your dataset( vector - format - vector - ensembl iD)
                            TestingSubsetNames, # "sig" module names
                            KEGGthres = 0.05, # significant level (default - 0.05)
                            species = "bta", 
                            id.type = "kegg",
                            Sig_list_out =Sig_list_out,
                            #biomart="ENSEMBL_MART_ENSEMBL",
                            #dataset="btaurus_gene_ensembl",
                            #host="http://www.ensembl.org",
                            #attributes = c("ensembl_gene_id","entrezgene_id"), # the items you need to retrive from the database
                            #filters="ensembl_gene_id", # with which keywords we match
                            keyword){
  #
  total_enrich = 0
  raw_pvalue_all = numeric()
  KEGG_results_b = list()
  KEGG_results_b_raw = list()
  library(biomaRt);library(gage);library(magrittr);library(stringr);library(ggplot2) # load functions source
  sdb = kegg.gsets(species = species, id.type = id.type, check.new = F) # get database1
  kegg.gs = sdb$kg.sets[sdb$sigmet.id] # organize database
  # get all genes in database
  # for splitting
  rexp <- "^(\\w+)\\s?(.*)$"
  genesKEGG = character()
  KEGGID = character()
  KEGGTERM = character()
  for (i in seq_along(names(kegg.gs))){
    tmp_gene = unlist(kegg.gs[i]);attributes(tmp_gene) = NULL
    KEGGID[i] = sub(rexp,"\\1",names(kegg.gs)[i])
    KEGGTERM[i] = sub(rexp,"\\2",names(kegg.gs)[i])
    genesKEGG = append(genesKEGG,tmp_gene,length(genesKEGG))
  }
  genesKEGG = unique(genesKEGG)
  # the output is a pdf and every single page will be the point plot of the enriched item of a specific module.
  # pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    if (i%%1==0){message("Now digging in module #",i)} # can change the 
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesKEGG])
    S = length(sig.genes[sig.genes %in% genesKEGG]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ID=character(),
                     Term=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    # Double loop: trying to go through every single KEGG.db, so extract each one first
    for (j in 1:length(names(kegg.gs))){
      #KEGG_Index = unlist(str_split(names(kegg.gs)[j]," ",2))[1] # split to get the GO-index
      #KEGG_Name = unlist(str_split(names(kegg.gs)[j]," ",2))[2] # split to get the Go name
      all_ENTER_temp = (as.vector(unlist(kegg.gs[j]))) #
      if (j%%100==0){message("checking on KEGG #",j,"-",KEGGID[j],"-",KEGGTERM[j])}
      # Calculate and overlap
      m = length(total.genes[total.genes %in% all_ENTER_temp]) # genes from target GO and in our dataset
      findG = sig.genes[sig.genes %in% all_ENTER_temp]
      s = length(findG) # # genes from target GO also in the non-preserved module
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 100)
      tmp = data.frame(ID= KEGGID[j], 
                       Term = KEGGTERM[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("KEGGID","KEGGTERM", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    KEGG_results_b_raw[[i]] = final_raw; names(KEGG_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched kegg raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= KEGGthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("KEGGID","KEGGTERM", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    KEGG_results_b[[i]] = final;names(KEGG_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched KEGG")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    # plotting
    #print(final %>% 
    #        top_n(dim(final)[1], wt= -pvalue) %>% 
    #        mutate(hitsPerc=Significant_Genes*100/Total_Genes)%>%
    #        ggplot(aes(x=hitsPerc,
    #                   y=KEGG_Name,
    #                   colour=pvalue,
    #                   size=Significant_Genes)) +
    #        xlim(0,max(final$hitsPerc)+5)+
    #        geom_point() +
    #        theme_gray()+
    #        labs(title= paste("KEGG Enrichment in module",module_name), x="Hits (%)", y="Kegg term", colour="p value", size="Count")+
    #       theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  }
  #dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  save(KEGG_results_b,
       KEGG_results_b_raw,
       raw_pvalue_distribution, 
       file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significantly pathways found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",KEGGthres)
  message("Nice! - KEGG enrichment finished and data saved")}
#########################################################################################################################
ConvertNformat = function(bg_gene,
                          TestingSubsetNames,
                          TestingModAssign,
                          keyword = "Ensembl2Entrez_Convert"){
  # Get match information
  key.symbol = AnnotationDbi::keys(org.Bt.eg.db,  keytype = c("ENSEMBL"))
  entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol), 
                                         columns = c("ENTREZID"),keytype = "ENSEMBL") %>% 
    dplyr::distinct(ENSEMBL,.keep_all= TRUE)
  #
  library(tidyverse)
  Gather_all = data.frame(ENSEMBL  =  bg_gene,
                          assign = TestingModAssign) %>% 
    dplyr::left_join(entrezUniverse, by  = c("ENSEMBL" = "ENSEMBL"))
  names( Gather_all)[3] = "ENTREZID_final"
  #
  Sig_list_out = list();Total_list_out = list()
  Sig_list_out_entrez = list();Total_list_out_entrez = list()
  Sig_list_out_ens = list();Total_list_out_ens = list()
  for (i in seq_along(TestingSubsetNames)){
    #
    target = TestingSubsetNames[i]
    tmp01 = dplyr::filter(Gather_all,assign == target)
    names(tmp01)[3] = "ENTREZID_final"
    Sig_list_out[[i]] = tmp01;names(Sig_list_out)[i] = TestingSubsetNames[i]
    tmp1 = dplyr::select(tmp01,ENTREZID_final) %>% dplyr::distinct() %>% na.omit();attributes(tmp1) = NULL
    Sig_list_out_entrez[[i]] = tmp1
    names(Sig_list_out_entrez)[i] = TestingSubsetNames[i]
    tmp2 = dplyr::select(tmp01,ENSEMBL) %>% dplyr::distinct() %>% na.omit();attributes(tmp1) = NULL
    Sig_list_out_ens[[i]] = tmp2;names(Sig_list_out_ens)[i] = TestingSubsetNames[i]
    names(Sig_list_out_entrez)[i] = TestingSubsetNames[i]
    
  }
  Total_list_out_tmp = list(unique(Gather_all$ENTREZID_final));names(Total_list_out_tmp) = "total_genes_entrez"
  Total_list_out_entrez = rep(Total_list_out_tmp, length(Sig_list_out_entrez))
  Total_list_out_tmp2 = list(unique(Gather_all$ENSEMBL));names(Total_list_out_tmp2) = "total_genes_ens"
  Total_list_out_ens = rep(Total_list_out_tmp2, length(Sig_list_out_ens))
  
  save(Sig_list_out,
       Sig_list_out_entrez,Total_list_out_entrez,
       Sig_list_out_ens,Total_list_out_ens,
       file = paste(trimws(keyword),".RData",sep = ""))
  message("Nice! Conversion finished")
}

print("update 0108 4pm")
