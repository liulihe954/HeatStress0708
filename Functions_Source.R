#===========================================================================================
#                                0. Package preparations                                  ##
#===========================================================================================
#devtools::install_github("jakesherman/easypackages")
library(easypackages)
my_packages <- c("sva","WGCNA","ppcor","edgeR","clusterProfiler","magrittr",
                 "dplyr", "ggplot2", "biomaRt","gage","doParallel",
                 "limma","recount","pamr","stringr")
libraries(my_packages)
#===========================================================================================
#                                1. "Massage" - Data processing                           ##
#===========================================================================================
# Classical way, networkData(rawdata,genes in rows and samples in cols), 
#                 n1-n2 (number of ref and treatmnent group)
#                 perct (the percentage to REMOVE based on VARIANCE across two conditions)            
DataPre   = function(networkData, cousin = 0.4, n1, n2, perct){
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
  #dim(networkData_normalized)
  # step 3 - log 2 trans
  # log trans
  #networkData_log2 = log2(networkData_normalized+2)
  # step 4 - filter out bottom xx% variation
  # select most var
  networkData_normalized$variance = apply(networkData_normalized, 1, var)
  networkData_50var = networkData_normalized[networkData_normalized$variance >= quantile(networkData_normalized$variance,c(perct)), ] #50% most variable genes
  networkData_50var$variance <- NULL
  #dim(networkData14_log2_50var)
  # step 5 - pca correction
  networkData_final = networkData_50var
  #save files for out put
  save(networkData_final,
       networkData_normalized_normfactors,
       networkData_normalized,
       file = paste(deparse(substitute(networkData)),"prepare with corrections","_top",100*(1-perct),".RData",sep = ""))
  return(list(Processed_final = networkData_final))
}
# Same rationales but FANCY way, remove confounding artifacts
# (ref - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1700-9 )
DataPre_C = function(networkData, cousin = 0.4, n1, n2, perct){
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
  #dim(networkData_normalized)
  # step 3 - log 2 trans
  # log trans
  networkData_log2 = log2(networkData_normalized+2)
  # step 4 - filter out bottom xx% variation
  # select most var
  networkData_log2$variance = apply(networkData_log2,1,var)
  networkData_log2_50var = networkData_log2[networkData_log2$variance >= quantile(networkData_log2$variance,c(perct)),] 
  networkData_log2_50var$variance <- NULL
  #dim(networkData14_log2_50var)
  # step 5 - pca correction 
  networkData_correction = Correct_pca(networkData_log2_50var,"leek")
  networkData_final = data.frame(networkData_correction$exprs_corrected_norm); 
  names(networkData_final) = names(networkData_log2_50var)
  save(networkData_final,
       networkData_log2_50var,
       networkData_normalized_normfactors,
       networkData_normalized,
       file = paste(deparse(substitute(networkData)),"prepare with corrections","_top",100*(1-perct),".RData",sep = ""))
  return(list(Corrected_log2_PC = networkData_final))
}
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
##########################################################################################
Kegg_Enrich_Plot = function(ENS_ID_all, # all genes in your dataset( vector - format - vector - ensembl iD)
                            TestingGroupAssignment, # color/module assignments 
                            TestingSubsetNames, # "sig" module names
                            KEGGthres = 0.05, # significant level (default - 0.05)
                            species = "bta", 
                            id.type = "kegg",
                            biomart="ENSEMBL_MART_ENSEMBL",
                            dataset="btaurus_gene_ensembl",
                            host="http://www.ensembl.org",
                            attributes = c("ensembl_gene_id","entrezgene_id"), # the items you need to retrive from the database
                            filters="ensembl_gene_id", # with which keywords we match
                            keyword){ # keyword is just for easy file naming, the keyword you provide will show as the main part as the output file name. e.g. "Day14_bicor_c_enrich"
  ##############           Matching (ens -> entrez); then plotting     #############################        
  # Load database from the downtown and try to match
  total_enrich = 0 # for counting the total enriched pathways
  library(biomaRt);library(gage);library(magrittr);library(stringr) # load functions source
  sdb = kegg.gsets(species = species, id.type = id.type, check.new = F) # get database1
  kegg.gs = sdb$kg.sets[sdb$sigmet.id] # organize database
  #length(sdb$kg.sets);str(kegg.gs)
  # Note, in kegg related database, the identifier is EntrezID, so we need to convert EnsemblIDs.
  mart <- biomaRt::useMart(biomart = biomart,
                           dataset = dataset,
                           host = host) # get database 
  # with annotation (matching in the databset), it actually the MATCHIG (intersection) with our dataset.
  annot_all <- getBM(attributes = attributes,
                     filters = filters,
                     values = ENS_ID_all,
                     mart = mart) 
  # extract EnsemblID for later use
  ENS_ID_all_annot <- as.vector(annot_all[,1])
  # extract EntrezID for later use
  ENTER_ID_all_annot <- as.vector(annot_all[,2])
  #length(ENTER_ID_all_annot);length(ENTER_ID_all_annot)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("We have ",length(ENS_ID_all)," Ensembl IDs ",
          " And we have ",length(ENS_ID_all_annot)," EntreZ IDs")
  if (length(ENS_ID_all)- length(ENS_ID_all_annot)>0){
    message("Note!! - Matching Done! we lost ",(-length(ENS_ID_all) + length(ENS_ID_all_annot))," nodes.")
  } else if (length(ENS_ID_all)- length(ENS_ID_all_annot) < 0){
    message("Note!! - Matching Done! we gained ",(-length(ENS_ID_all) + length(ENS_ID_all_annot))," nodes.")
  } else (message(" Matching Done! Perfectly matched! "))
   ###############################Matching done; now plotting #################################
  # the output is a pdf and every single page will be the point plot of the enriched item of a specific module.
  pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    if (i%%5==0){message("Now digging in module #",i)} # can change the 
    module_name = TestingSubsetNames[i]# get non preserved modules / (nodes)
    nopresID = as.vector(ENS_ID_all[which(TestingGroupAssignment == module_name)]) # Matching every module
    nopresENTER = annot_all[ENS_ID_all_annot %in% nopresID,2] # Note. multiple EntrezID could point to single EnsemblID!!
    N = length(ENTER_ID_all_annot) # Big N - overlap of our dataset(top20) and the KEGG.db
    S = length(ENS_ID_all_annot[ENS_ID_all_annot%in%nopresID]) # Big S, find those in both top 20 and no-pres modules
    out = data.frame(KEGG=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric()) # formating
    # Double loop: trying to go through every single KEGG.db, so extract each one first
    for (j in 1:length(names(kegg.gs))){
      KEGG_Index = unlist(str_split(names(kegg.gs)[j]," ",2))[1] # split to get the GO-index
      KEGG_Name = unlist(str_split(names(kegg.gs)[j]," ",2))[2] # split to get the Go name
      all_ENTER_temp = (as.vector(unlist(kegg.gs[j]))) #
      if (j%%20==0){message("checking on KEGG #",j,"-",KEGG_Index,"-",KEGG_Name)}
      # Calculate and overlap
      m = length(ENTER_ID_all_annot[ENTER_ID_all_annot %in% all_ENTER_temp]) # genes from target GO and in our dataset
      s = length(nopresENTER[nopresENTER %in% all_ENTER_temp]) # # genes from target GO also in the non-preserved module
      # format a matrix
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
      tmp = data.frame(KEGG = KEGG_Index, Name = KEGG_Name, totalG = m, sigG = s, Pvalue = Pval)
      out = rbind(out,tmp)}
    # select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue < KEGGthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("KEGG.ID","KEGG_Name", "Total_Genes", "Significant_Genes", "pvalue")
    final = final %>% top_n(dim(final)[1], wt= -pvalue) %>% mutate(hitsPerc=Significant_Genes*100/Total_Genes)
    KEGG_results_b[[i]] = final
    total_enrich = total_enrich + nrow(final)
    # plotting
    print(final %>% 
          top_n(dim(final)[1], wt= -pvalue) %>% 
          mutate(hitsPerc=Significant_Genes*100/Total_Genes)%>%
          ggplot(aes(x=hitsPerc,
                       y=KEGG_Name,
                       colour=pvalue,
                       size=Significant_Genes)) +
          xlim(0,max(final$hitsPerc)+5)+
          geom_point() +
          theme_gray()+
          labs(title= paste("KEGG Enrichment in module",module_name), x="Hits (%)", y="Kegg term", colour="p value", size="Count")+
          theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
          theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  }
  dev.off()
  save(KEGG_results_b, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significantly pathways found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",KEGGthres)
  message("Nice! - KEGG enrichment finished and data saved")
  Match_Annot=list(All_Ens_ID_after_annot=ENS_ID_all_annot,
                   All_Entrez_ID_after_annot=ENTER_ID_all_annot,
                   Enrich_Compile = KEGG_results_b)
  return(Match_Annot)}

#####################################################################################
#=== Following are for testing 'wgcna' --- ignore unless you find them relevent ===#
#####################################################################################
#ENS_ID_all <- colnames(datExpr14_cl)
#nonpres_index_b = (which(Zsummary14_b < 2))
#nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
#nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
#nonpres_modulenames_b = nonpres_modulenames_b[c(15,17)]
#KEGG_results_b = list()
#TestingKegg = Kegg_Enrich_Plot(ENS_ID_all,
#                               KEGGthres = 0.05,
#                               TestingGroupAssignment = moduleColors14_b_cl, 
#                               TestingSubsetNames = nonpres_modulenames_b,
#                               keyword = "KEGG_Enrichment_just_testing")
#####################################################################################




#===========================================================================================
#                              3.Gene Ontology enrichment                                   ##
#===========================================================================================
##########################################################################################
### Function necessities: Same rationales with KEGG enrich 
###                                    ---
###                        !but here we dont need convert from EnsemblID to EntrezID!
##########################################################################################
Go_Enrich_Plot = function(total.genes = total.genes,
                          TestingGroupAssignment, 
                          TestingSubsetNames,
                          GOthres = 0.05,
                          biomart="ensembl",
                          dataset="btaurus_gene_ensembl",
                          host="http://www.ensembl.org",
                          attributes = c("ensembl_gene_id", "external_gene_name","go_id","name_1006"),
                          keyword = "GO_Enrichment_in_modules_bicor_c_day14_new_testing"){
  # 
  total_enrich = 0
  library(biomaRt);library(gage);library(magrittr)# load pkg
  ## Analysis bosTau annotation: GO
  database2 <- biomaRt::useMart(biomart=biomart,
                                dataset=dataset,
                                host=host)
  gene2 <- getBM(attributes = attributes,
                 mart = database2)
  #dim(gene2); length(unique(gene2$ensembl_gene_id)); length(unique(gene2$go_id))
  goName = unique(gene2[,c(3,4)]);goName = goName[order(goName$go_id),];goName = goName[-1,]
  GO = goName$go_id;Name = goName$name_1006
  genesGO = unique(subset(gene2,go_id != "")$ensembl_gene_id)
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of GO sets to check: ",length(GO)," with total number of names: ",length(Name))
  # plot
  pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on module #",i," - ",TestingSubsetNames[i])
    module_name = TestingSubsetNames[i]
    nopresID_GO = as.vector(colnames(datExpr14_cl)[which(TestingGroupAssignment == module_name)])
    sig.genes = nopresID_GO # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesGO])
    S = length(sig.genes[sig.genes %in% genesGO]) #
    out = data.frame(GO=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric())
    message("Module size of ",TestingSubsetNames[i],": ", length(nopresID_GO))
    for(j in 1:length(GO)){
      if (j%%100 == 0) {message("tryingd on GO ",j," - ",GO[j]," - ",Name[j])}
      gENEs = subset(gene2, go_id == GO[i])$ensembl_gene_id # all gene in target GO
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target GO also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
      tmp = data.frame(GO = GO[j], Name = Name[j], totalG = m, sigG = s, Pvalue = Pval)
      out = rbind(out,tmp)}
    # select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue < GOthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue")
    final = final %>% top_n(dim(final)[1], wt= -pvalue)%>%mutate(hitsPerc = Significant_Genes*100/Total_Genes)
    GO_results_b[[i]] = final;names(GO_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched GO")
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    #
    print(final %>%
          top_n(dim(final)[1], wt= -pvalue)%>%
          mutate(hitsPerc = Significant_Genes*100/Total_Genes) %>% ## signi genes, v1 = all genes in the go.
          ggplot(aes(x = hitsPerc,
                       y = GO_Name,
                       colour = pvalue,
                       size = Significant_Genes)) +
            #xlim(0,)+
          geom_point() +
          theme_gray()+
          labs(title= paste("GO Enrichment in module",module_name), x="Hits (%)", y="GO term", colour="p value", size="Count")+
          theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
          theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
          theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  }
  dev.off()
  message(total_enrich," significant GO terms found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",GOthres)
  message("Nice! - GO enrichment finished and data saved")
  save(GO_results_b, file = paste(trimws(keyword),".RData",sep = ""))}

#####################################################################################
#=== Following are for testing 'wgcna' --- ignore unless you find them relevent ===#
#####################################################################################
#total.genes = colnames(datExpr14_cl) # total genes in your dataset
#nonpres_index_b = (which(Zsummary14_b < 2))
#nonpres_modulenames_b = rownames(Z.PreservationStats14_b)[nonpres_index_b]
#nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)];GO_results_b = list() 
#nonpres_modulenames_b = nonpres_modulenames_b[c(13,19)]
#total.genes = colnames(datExpr14_cl) # total genes in your dataset
#TestingGO = Go_Enrich_Plot(total.genes,
#                           GOthres = 0.05,
#                           TestingGroupAssignment = moduleColors14_b_cl, 
#                           TestingSubsetNames = nonpres_modulenames_b,
#                           keyword = "GO_Enrichment_just_testing")
#####################################################################################