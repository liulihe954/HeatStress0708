require(sva);require(WGCNA)
require(ppcor);require(dplyr)
require(edgeR);require(clusterProfiler)
require(ggplot2);require(magrittr)
require(biomaRt);require(gage);require(doParallel)
require(limma);require(recount);require(pamr)
#================================================================================================
###                                       1. dataprep                                      ######
#================================================================================================
## read
#setwd("/Users/liulihe95/Desktop/HeatStress0708");getwd()
options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
CowsID_ht = c("6334","8514","8971","8867","8841","8966")# ID of heat group
CowsID_cl = c("8252","8832","8896","8983","8897","8862")# ID of cool group
# (all the IDs) 27570 genes; 12 cows; 3 time points 
networkData = read.csv("count_matrix.tsv",sep = "\t",header = T);rownames(networkData) = networkData[,1];networkData = networkData[,-1] 
dim(networkData)
## find index 
# day 14
names14_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos14 = grep(".14",names14_last3)
names14 = names(networkData)[pos14];pos14_ht = substr(names14,2,5) %in% CowsID_ht;pos14_cl = substr(names14,2,5) %in% CowsID_cl
# day 42
names42_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos42 = grep(".42",names42_last3)
names42 = names(networkData)[pos42];pos42_ht = substr(names42,2,5) %in% CowsID_ht;pos42_cl = substr(names42,2,5) %in% CowsID_cl
# day 84
names84_last3 = substr(names(networkData), nchar(names(networkData)[2])-2, nchar(names(networkData)[2])); pos84 = grep(".84",names84_last3)
names84 = names(networkData)[pos84];pos84_ht = substr(names84,2,5) %in% CowsID_ht;pos84_cl = substr(names84,2,5) %in% CowsID_cl
#compose index
column_14_cl = pos14[pos14_cl];column_14_ht = pos14[pos14_ht]
column_42_cl = pos42[pos42_cl];column_42_ht = pos42[pos42_ht]
column_84_cl = pos84[pos84_cl];column_84_ht = pos84[pos84_ht]
# -4 in day14; -4 in day 42; -7 in day84 (think about which one to respect)
networkData14 = networkData[,c(column_14_cl,column_14_ht)]
networkData42 = networkData[,c(column_42_cl,column_42_ht)]
networkData84 = networkData[,c(column_84_cl,column_84_ht)]
dim(networkData14);dim(networkData42);dim(networkData84)

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



#=======================================day14=========================================================
networkData14_perc20 = DataPre(networkData14, cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData14_perc20_C = DataPre_C(networkData14,cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData14_final <- networkData14_perc20$Processed_final
networkData14_final_c <- networkData14_perc20_C$Corrected_log2_PC

datExpr14_cl = t(networkData14_final[,colnames(networkData14_final) %in% names(networkData)[column_14_cl] ])
datExpr14_ht = t(networkData14_final[,colnames(networkData14_final) %in% names(networkData)[column_14_ht] ])
datExpr14_cl = data.frame(datExpr14_cl);datExpr14_ht = data.frame(datExpr14_ht)

datExpr14_cl_c = t(networkData14_final_c[,colnames(networkData14_final_c) %in% names(networkData)[column_14_cl] ])
datExpr14_ht_c = t(networkData14_final_c[,colnames(networkData14_final_c) %in% names(networkData)[column_14_ht] ])
datExpr14_cl_c = data.frame(datExpr14_cl_c);datExpr14_ht_c = data.frame(datExpr14_ht_c)

#=======================================day42=========================================================
networkData42_perc20 = DataPre(networkData42, cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData42_perc20_C = DataPre_C(networkData42,cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData42_final <- networkData42_perc20$Processed_final
networkData42_final_c <- networkData42_perc20_C$Corrected_log2_PC

datExpr42_cl = t(networkData42_final[,colnames(networkData42_final) %in% names(networkData)[column_42_cl] ])
datExpr42_ht = t(networkData42_final[,colnames(networkData42_final) %in% names(networkData)[column_42_ht] ])
datExpr42_cl = data.frame(datExpr42_cl);datExpr42_ht = data.frame(datExpr42_ht)

datExpr42_cl_c = t(networkData42_final_c[,colnames(networkData42_final_c) %in% names(networkData)[column_42_cl] ])
datExpr42_ht_c = t(networkData42_final_c[,colnames(networkData42_final_c) %in% names(networkData)[column_42_ht] ])
datExpr42_cl_c = data.frame(datExpr42_cl_c);datExpr42_ht_c = data.frame(datExpr42_ht_c)

#=======================================day84=========================================================
networkData84_perc20 = DataPre(networkData84, cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData84_perc20_C = DataPre_C(networkData84,cousin = 0.4, n1 = 6, n2 = 6, perct = .80)
networkData84_final <- networkData84_perc20$Processed_final
networkData84_final_c <- networkData84_perc20_C$Corrected_log2_PC

datExpr84_cl = t(networkData84_final[,colnames(networkData84_final) %in% names(networkData)[column_84_cl] ])
datExpr84_ht = t(networkData84_final[,colnames(networkData84_final) %in% names(networkData)[column_84_ht] ])
datExpr84_cl = data.frame(datExpr84_cl);datExpr84_ht = data.frame(datExpr84_ht)

datExpr84_cl_c = t(networkData84_final_c[,colnames(networkData84_final_c) %in% names(networkData)[column_84_cl] ])
datExpr84_ht_c = t(networkData84_final_c[,colnames(networkData84_final_c) %in% names(networkData)[column_84_ht] ])
datExpr84_cl_c = data.frame(datExpr84_cl_c);datExpr84_ht_c = data.frame(datExpr84_ht_c)

#================================================================================================

powers = c(c(1:10), seq(from = 12, to=30, by=2))
#================================================================================================
###                                     2. day14 cor                                       ######
#================================================================================================
sft_cl_14_cor = pickSoftThreshold(datExpr14_cl, powerVector = powers, verbose = 0)
#================================================================================================
###                                     3. day14 cor_c                                     ######
#================================================================================================
sft_cl_14_cor_c = pickSoftThreshold(datExpr14_cl_c, powerVector = powers, verbose = 0)
#================================================================================================
###                                     4. day14 bicor                                     ######
#================================================================================================
sft_cl_14_bicor = pickSoftThreshold(datExpr14_cl, powerVector = powers, corFnc = "bicor",verbose = 0)
#================================================================================================
###                                    5. day14 bicor_c                                    ######
#================================================================================================
sft_cl_14_bicor_c = pickSoftThreshold(datExpr14_cl_c, powerVector = powers, corFnc = "bicor",verbose = 0)

save(sft_cl_14_cor,
     sft_cl_14_cor_c,
     sft_cl_14_bicor,
     sft_cl_14_bicor_c,
     file = "soft_thres_four_comb_day14.RData")


#================================================================================================
###                                     2. day42 cor                                       ######
#================================================================================================
sft_cl_42_cor = pickSoftThreshold(datExpr42_cl, powerVector = powers, verbose = 0)
#================================================================================================
###                                     3. day14 cor_c                                     ######
#================================================================================================
sft_cl_42_cor_c = pickSoftThreshold(datExpr42_cl_c, powerVector = powers, verbose = 0)
#================================================================================================
###                                     4. day14 bicor                                     ######
#================================================================================================
sft_cl_42_bicor = pickSoftThreshold(datExpr42_cl, powerVector = powers, corFnc = "bicor",verbose = 0)
#================================================================================================
###                                    5. day14 bicor_c                                    ######
#================================================================================================
sft_cl_42_bicor_c = pickSoftThreshold(datExpr42_cl_c, powerVector = powers, corFnc = "bicor",verbose = 0)
#================================================================================================
save(sft_cl_42_cor,
     sft_cl_42_cor_c,
     sft_cl_42_bicor,
     sft_cl_42_bicor_c,
     file = "soft_thres_four_comb_day42.RData")

#================================================================================================
###                                     2. day84 cor                                       ######
#================================================================================================
sft_cl_84_cor = pickSoftThreshold(datExpr84_cl, powerVector = powers, verbose = 0)
#================================================================================================
###                                     3. day14 cor_c                                     ######
#================================================================================================
sft_cl_84_cor_c = pickSoftThreshold(datExpr84_cl_c, powerVector = powers, verbose = 0)
#================================================================================================
###                                     4. day14 bicor                                     ######
#================================================================================================
sft_cl_84_bicor = pickSoftThreshold(datExpr84_cl, powerVector = powers, corFnc = "bicor",verbose = 0)
#================================================================================================
###                                    5. day14 bicor_c                                    ######
#================================================================================================
sft_cl_84_bicor_c = pickSoftThreshold(datExpr84_cl_c, powerVector = powers, corFnc = "bicor",verbose = 0)
#================================================================================================
save(sft_cl_84_cor,
     sft_cl_84_cor_c,
     sft_cl_84_bicor,
     sft_cl_84_bicor_c,
     file = "soft_thres_four_comb_day84.RData")

