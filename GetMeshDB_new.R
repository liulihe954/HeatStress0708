# Save mesh for later use
library(MeSH.db);library(MeSH.Bta.eg.db)
library(tidyverse);library(gage);library(magrittr)
#
#keyword_outer = "MeshDB"
#DB = paste(keyword_outer,".RData",sep = "")
#load(DB)
#
#library(ggplot2);library(biomaRt) # load pkg
# raw data for retrive MESHid and all details linked
#KEY = keys(MeSH.db, keytype = "MESHID")
#List = MeSHDbi::select(MeSH.db, keys = KEY[5], columns = columns(MeSH.db), keytype = "MESHID")
#List = select(MeSH.db, keys = KEY[1:3], columns = columns(MeSH.db), keytype = "MESHID")
#Match_List = dplyr::select(List, MESHID, MESHTERM)

#library(MeSH.Bta.eg.db)
#key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
#list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
#  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% c("D","G")) %>% 
#  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
#save(KEY,List,Match_List,key_Bta,list_Bta,file = "MeshDB.RData")
load("MeshDB.RData")
##
library(tidyverse)
MeshRecords = unique(list_Bta[,c("MESHID","MESHTERM")]) %>% arrange(MESHID)
MeshID = na.omit(MeshRecords$MESHID)
MeshTerm = na.omit(MeshRecords$MESHTERM)

DB_List = list()
for ( p in seq_along(MeshID)){
  tmp = subset(list_Bta, MESHID == MeshID[p])$GENEID
  DB_List[[p]] = tmp #
  names(DB_List)[p]  <- paste(MeshID[p],"-",MeshTerm[p])
}

save(DB_List,KEY,List,Match_List,key_Bta,list_Bta,file = "MeshDB_new.RData")

