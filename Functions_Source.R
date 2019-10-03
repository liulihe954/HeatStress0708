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
### Function output: a. a plot with all the enriched items, one page for each module.
###                  b. .RData containing a long list, each element show all the enriched items for corresponding module.
###                      same length with that of all non-preserved module
###
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
                            keyword){
  #
  total_enrich = 0
  raw_pvalue_all = numeric()
  KEGG_results_b = list()
  KEGG_results_b_raw = list()
  library(biomaRt);library(gage);library(magrittr);library(stringr);library(ggplot2) # load functions source
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
                     mart = mart) %>% tidyr::drop_na(.)
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
    #module_name = TestingSubsetNames[i]# get non preserved modules / (nodes)
    nopresID = as.vector(ENS_ID_all[which( TestingGroupAssignment == module_name)]) # Matching every module
    nopresENTER = annot_all[ENS_ID_all_annot %in% nopresID,2] # Note. multiple EntrezID could point to single EnsemblID!!
    N = length(ENTER_ID_all_annot) # Big N - overlap of our dataset(top20) and the KEGG.db
    S = length(ENTER_ID_all_annot[ENTER_ID_all_annot%in%nopresENTER]) # Big S, find those in both top 20 and no-pres modules
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
    # put the pvalues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # select those has 4 more gene in common and pvalue smaller than 0.05
    # selection vs raw
    ot = subset(out,totalG > 4 & Pvalue < KEGGthres) # select
    ot_raw = out # no-select; keep all
    # 
    final = ot[order(ot$Pvalue),];colnames(final) = c("KEGG.ID","KEGG_Name", "Total_Genes", "Significant_Genes", "pvalue")
    final = final %>% top_n(dim(final)[1], wt= -pvalue) %>% mutate(hitsPerc=Significant_Genes*100/Total_Genes)
    KEGG_results_b[[i]] = final
    total_enrich = total_enrich + nrow(final)
    # 
    final_raw = ot_raw[order(ot_raw$Pvalue),];colnames(final_raw) = c("KEGG.ID","KEGG_Name", "Total_Genes", "Significant_Genes", "pvalue")
    final_raw = final_raw %>% top_n(dim(final_raw)[1], wt= -pvalue) %>% mutate(hitsPerc=Significant_Genes*100/Total_Genes)
    KEGG_results_b_raw[[i]] = final_raw
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
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){
    raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))
  }
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  save(ENS_ID_all_annot,
       ENTER_ID_all_annot,
       KEGG_results_b,
       KEGG_results_b_raw,
       raw_pvalue_all, 
       file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significantly pathways found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",KEGGthres)
  message("Nice! - KEGG enrichment finished and data saved")
  Match_Annot=list(All_Ens_ID_after_annot=ENS_ID_all_annot,
                   All_Entrez_ID_after_annot=ENTER_ID_all_annot,
                   Enrich_Compile = KEGG_results_b,
                   Enrich_Compile_raw = KEGG_results_b_raw,
                   All_pvalue_distribution =  raw_pvalue_distribution)
  return(Match_Annot)}
#####################################################################################
Parse_KEGG_Results = function(KEGG_results_b){
  all_enrich_KEGG = data.frame(ID=character(),
                               Description=character(),
                               GeneRatio=character(),
                               BgRatio=character(),
                               pvalue=numeric(),
                               p.adjust=numeric(),
                               qvalue=numeric(),
                               geneID=character(),
                               Count=numeric(),
                               stringsAsFactors=FALSE)
  for (i in 1:length(KEGG_results_b)){
    len = dim(data.frame(KEGG_results_b[i]))[1]
    if (len> 0){
      all_enrich_KEGG = rbind(all_enrich_KEGG,data.frame(KEGG_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_KEGG)[1]
  total_modules = length(KEGG_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_KEGG)
}
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
                          attributes = c("external_gene_name","go_id","name_1006"),
                          keyword = "GO_Enrichment_in_modules_bicor_c_day14_new_testing"){
  # 
  total_enrich = 0
  raw_pvalue_all = numeric()
  GO_results_b = list()
  GO_results_b_raw = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr)# load pkg
  ## Analysis bosTau annotation: GO
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  goName = unique(gene[,c(2,3)]);goName = goName[order(goName$go_id),];goName = goName[-1,]
  GO = goName$go_id
  Name = goName$name_1006
  #length(GO)
  genesGO = unique(subset(gene,go_id != "")$ensembl_gene_id)[-1]
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of GO sets to check: ",length(GO)," with total number of names: ",length(Name))
  # plot
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on module #",i," - ",TestingSubsetNames[i])
    module_name = TestingSubsetNames[i]
    nopresID_GO = as.vector(colnames(datExpr14_cl)[which(TestingGroupAssignment == module_name)])
    sig.genes = nopresID_GO # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesGO])
    S = length(sig.genes[sig.genes %in% genesGO]) #
    ExternalLoss_total = paste((length(total.genes) - N),(length(total.genes) - N)/N,sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),(length(sig.genes) - S)/S,sep = "/")
    out = data.frame(GO=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(nopresID_GO))
    for(j in 1:length(GO)){
      if (j%%100 == 0) {message("tryingd on GO ",j," - ",GO[j]," - ",Name[j])}
      gENEs = subset(gene, go_id == GO[i])$ensembl_gene_id # all gene in target GO
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target GO also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
      tmp = data.frame(GO = GO[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    #ot_raw = subset(out,totalG > 4 & Pvalue < GOthres)
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% top_n(dim(final_raw)[1], wt= -pvalue_r)%>%mutate(hitsPerc = Significant_Genes*100/Total_Genes)
    GO_results_b_raw[[i]] = final_raw;names(GO_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched GO raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue < GOthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue","ExternalLoss_total","InternalLoss_sig")
    final = final %>% top_n(dim(final)[1], wt= -pvalue)%>%mutate(hitsPerc = Significant_Genes*100/Total_Genes)
    GO_results_b[[i]] = final;names(GO_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched GO")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    #
    #print(final %>%
    #        top_n(dim(final)[1], wt= -pvalue)%>%
    #        mutate(hitsPerc = Significant_Genes*100/Total_Genes) %>% ## signi genes, v1 = all genes in the go.
    #        ggplot(aes(x = hitsPerc,
    #                   y = GO_Name,
    #                   colour = pvalue,
    #                   size = Significant_Genes)) +
            #xlim(0,)+
    #        geom_point() +
    #        theme_gray()+
    #        labs(title= paste("GO Enrichment in module",module_name), x="Hits (%)", y="GO term", colour="p value", size="Count")+
    #        theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  }
  #dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){
    raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))
  }
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  raw_pvalue_distribution
  save(GO_results_b, GO_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant GO terms found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",GOthres)
  message("Nice! - GO enrichment finished and data saved")}


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
#######################################################################################
#   Function Parse_GO funtion, --- unlist the list and put them in one data.frame   #
#######################################################################################
Parse_GO_Results = function(GO_results_b){
  all_enrich_GO = data.frame(ID=character(),
                             Description=character(),
                             Total_gene=numeric(),
                             Significant_gene=numeric(),
                             pvalue=numeric(),
                             ExternalLoss_total = character(),
                             InternalLoss_total = character(),
                             HitPerc = numeric(),
                             stringsAsFactors=FALSE)
  for (i in 1:length(GO_results_b)){
    len = dim(data.frame(GO_results_b[i]))[1]
    if (len> 0){
      tmp = data.frame(GO_results_b[i])
      names(tmp) = names(all_enrich_GO)
      all_enrich_GO = rbind(all_enrich_GO,tmp)
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_GO)[1]
  total_modules = length(GO_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_GO)
}

#######################################################################################
#         Function Semantic similarity plotting, calculate correlation score          #
#######################################################################################
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
#########################################################################################################################
Parse_Results = function(KEGG_results_b){
  all_enrich_KEGG = data.frame(ID=character(),
                               Description=character(),
                               GeneRatio=character(),
                               BgRatio=character(),
                               pvalue=numeric(),
                               p.adjust=numeric(),
                               qvalue=numeric(),
                               geneID=character(),
                               Count=numeric(),
                               stringsAsFactors=FALSE)
  for (i in 1:length(KEGG_results_b)){
    len = dim(data.frame(KEGG_results_b[i]))[1]
    if (len> 0){
      all_enrich_KEGG = rbind(all_enrich_KEGG,data.frame(KEGG_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_KEGG)[1]
  total_modules = length(KEGG_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_KEGG)
}

#########################################################################################################################
#########################################################################################################################
InterPro_Enrich = function(total.genes,
                           sig.genes,
                           TestingSubsetNames,
                           IPthres = 0.05,
                           biomart="ensembl",
                           dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           Identifier = "external_gene_name",
                           attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Interpro_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Interpro_results_b = list()
  Interpro_results_b_raw = list()
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
  if (Identifier == "ensembl_gene_id"){genesInterpro = unique(subset(gene,interpro != "")$ensembl_gene_id)
  } else if (Identifier == "external_gene_name") {
    genesInterpro= unique(subset(gene,interpro != "")$external_gene_name);
    genesInterpro = genesInterpro[-1]
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
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(Interpro)){
      if (j%%100 == 0) {message("tryingd on Interpro ",j," - ",Interpro[j]," - ",Name[j])}
      if (Identifier == "ensembl_gene_id"){
        gENEs = subset(gene, interpro == Interpro[j])$ensembl_gene_id
      } else if (Identifier == "external_gene_name") {
        gENEs = subset(gene, interpro == Interpro[j])$external_gene_name
      }
      m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target interpro also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Interpro = Interpro[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Interpro_results_b_raw[[i]] = final_raw; names(Interpro_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Interpro raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= IPthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
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
Parse_Interpro_Results = function(Interpro_results_b){
  all_enrich_Interpro = data.frame()
  for (i in 1:length(Interpro_results_b)){
    len = dim(data.frame(Interpro_results_b[i]))[1]
    if (len> 0){
      all_enrich_Interpro = rbind(all_enrich_Interpro,data.frame(Interpro_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_Interpro)[1]
  total_modules = length(Interpro_results_b)
  print(paste(total_hits,"hits found in",total_modules,"tested modules"))
  return(ParseResults = all_enrich_Interpro)
}
#########################################################################################################################
#########################################################################################################################
MESH_Enrich = function(total_genes_all,
                       sig_genes_all,
                       TestingSubsetNames,
                       Meshthres = 0.05,
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
  library(MeSH.db);library(MeSH.Bta.eg.db);library(tidyverse);library(gage);library(magrittr)
  library(ggplot2);library(biomaRt) # load pkg
  
  ### raw data for retrive MESHid and all details linked
  #
  #KEY = keys(MeSH.db, keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  ##List = select(MeSH.db, keys = KEY[1:3], columns = columns(MeSH.db), keytype = "MESHID")
  #Match_List = dplyr::select(List, MESHID, MESHTERM)
  ##head(Match_List) 
  ### Prepare Bta database
  #
  library(MeSH.Bta.eg.db)
  #key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
  #  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY == MeshCate) %>% 
  #  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
  #
  # alternatively
  keyword_outer = "MeshDB"
  DB = paste(keyword_outer,".RData",sep = "")
  load(DB)
  #
  # Get index
  genesMesh = unique(list_Bta$GENEID)
  MeshID = unique(list_Bta$MESHID)
  #MeshID = MeshID[1:1000]
  MeshTerm = unique(list_Bta$MESHTERM)
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Mesh to check: ",length(MeshID)," with total number of names: ",length(MeshTerm))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    # i = 2
    # sig.genes = Sig_list_out_entrez_test
    # total.genes = Total_list_out_entrez_test
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
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(MeshID)){
      if (j%%1000 == 0) {message("tryingd on MeshID ",j," - ",MeshID[j]," - ",MeshTerm[j])}
      #head(list_Bta)
      gENEs = subset(list_Bta, MESHID == MESHID[j])$GENEID
      m = length(total.genes[total.genes %in% gENEs]) # genes from target  and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target  also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(MeshID= MeshID[j], 
                       MeshTerm = MeshTerm[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Mesh_results_b_raw[[i]] = final_raw; names(Mesh_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Mesh raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Meshthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
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
Parse_Mesh_Results = function(Mesh_results_b){
  all_enrich_Mesh = data.frame()
  for (i in 1:length(Mesh_results_b)){
    len = dim(data.frame(Mesh_results_b[i]))[1]
    if (len> 0){
      all_enrich_Mesh = rbind(all_enrich_Mesh,data.frame(Mesh_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_Mesh)[1]
  total_modules = length(Mesh_results_b)
  print(paste(total_hits,"hits found in",total_modules,"tested modules"))
  return(ParseResults = all_enrich_Mesh)
}


#########################################################################################################################
#########################################################################################################################
Reactome_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           InputSource,
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
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg
  Reactome_gene =   unique(InputSource[,1])
  ReactomeID =      unique(InputSource[,2])
  ReactomeName =    unique(InputSource[,3])
  #
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Reactome to check: ",length(ReactomeID)," with total number of names: ",length(ReactomeName))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
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
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(ReactomeID)){
      if (j%%100 == 0) {message("tryingd on Reactome ",j," - ",ReactomeID[j]," - ",ReactomeName[j])}
      gENEs = subset(InputSource, ReactomeID == ReactomeID[j])$EntrezID
      m = length(total.genes[total.genes %in% gENEs]) 
      s = length(sig.genes[sig.genes %in% gENEs]) # 
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(ReactomeID = ReactomeID[j], 
                       ReactomeName = ReactomeName[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Reactome_results_b_raw[[i]] = final_raw; names(Reactome_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Reactomeid raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Reacthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
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
          " at the significance level of ",IPthres)
  message("Nice! - Reactome enrichment finished and data saved")}

#########################################################################################################################
Parse_Reactome_Results = function(Mesh_results_b){
  all_enrich_Mesh = data.frame()
  for (i in 1:length(Mesh_results_b)){
    len = dim(data.frame(Mesh_results_b[i]))[1]
    if (len> 0){
      all_enrich_Mesh = rbind(all_enrich_Mesh,data.frame(Mesh_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_Mesh)[1]
  total_modules = length(Mesh_results_b)
  print(paste(total_hits,"hits found in",total_modules,"tested modules"))
  return(ParseResults = all_enrich_Mesh)
} # did change name because its all the same
