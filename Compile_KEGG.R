load("KEGG_results_bicor_day14.RData")
load("KEGG_results_bicor_c_top20_day42.RData")
load("KEGG_results_bicor_c_day84.RData")
#
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
  total_hits = dim(all_enrich_KEGG)[1]
  total_modules = length(KEGG_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_KEGG)
}
KEGG_Results_Day14 = Parse_KEGG_Results(KEGG_results_b)
#
write.csv(KEGG_Results_Day14,"KEGG_Results_Day14.csv")
write.csv(KEGG_Results_Day14,"KEGG_Results_Day42.csv")
write.csv(KEGG_Results_Day14,"KEGG_Results_Day84.csv")

# day14   87 hits found in 59 non-preserved modules"
# day 42  20 hits found in 30 non-preserved modules"
# day 84  38 hits found in 38 non-preserved modules

#
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20")
#
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20/Day14")
KEGG_Results_Day14=read.csv("KEGG_Results_Day14.csv")
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20/Day42")
KEGG_Results_Day42=read.csv("KEGG_Results_Day42.csv")
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20/Day84")
KEGG_Results_Day84=read.csv("KEGG_Results_Day84.csv")
#
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20")

#install.packages("reshape")
library(reshape);library(reshape2)
list_allpath <- list(KEGG_Results_Day14[,c(1,3)],
                     KEGG_Results_Day42[,c(1,3)],
                     KEGG_Results_Day84[,c(1,3)])
merged_allpath <- merge_all(list_allpath)
#
AllPaths_index = as.vector(merged_allpath$X)
Compile_enrich_KEGG_info = cbind(Descriptions=as.vector(merged_allpath[,2]),
                                 Day14=rep(0,length(AllPaths_index)),
                                 Day42=rep(0,length(AllPaths_index)),
                                 Day84=rep(0,length(AllPaths_index)))
rownames(Compile_enrich_KEGG_info) = AllPaths_index


for (i in 1:length(AllPaths_index)){
  # for day14
  loci14 = which(AllPaths_index  %in% as.vector(KEGG_Results_Day14[,1]))
  Compile_enrich_KEGG_info[loci14,2] = "1"
  # for day42
  loci42 = which(AllPaths_index %in% as.vector(KEGG_Results_Day42[,1]))
  Compile_enrich_KEGG_info[loci42,3] = "1"
  # for day84
  loci84 = which(AllPaths_index  %in% as.vector(KEGG_Results_Day84[,1]))
  Compile_enrich_KEGG_info[loci84,4] = "1"
}

write.csv(Compile_enrich_KEGG_info,"Compile_enrich_KEGG_info.csv")
