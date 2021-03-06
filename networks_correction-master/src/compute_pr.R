rm(list = ls())

setwd("/work-zfs/abattle4/parsana/networks_correction/")

library(reshape2)
library(dplyr)
library(ggplot2)
library(recount)
library(gtools)

inputargs <- commandArgs(TRUE)
tiss <- inputargs[1]
pathways.fn <- inputargs[2]
plot.dir <- inputargs[3]
res.dir <- inputargs [4]
save.fn <- inputargs[5]
net_type <- inputargs[6]
path_shared <- eval(parse(text = inputargs[7]))

type.exp <- c("raw", "rin", "gc", "exonicRate", "mc", "quarterpc", "halfpc", "pc")

  # get gene ids - symbol mapping from recount dataset

  load("data/raw_subset.Rdata")
  dat.gene.symbol <- rowData(dat.expr[[tiss]])
  head(dat.gene.symbol)
  dat.gene.symbol <- sapply(dat.gene.symbol$symbol, function(x) x[1])
  dat.gene.symbol <- data.frame(gene_id = names(dat.gene.symbol), gene_symbol = dat.gene.symbol, stringsAsFactors = F)
  rm(dat.expr)

if(path_shared){
    true.positive.list <- read.delim(pathways.fn, header = T, row.names = 1, stringsAsFactors = F, sep = " ")
    true.positive.list <- true.positive.list[true.positive.list[,1] %in% dat.gene.symbol$gene_symbol & true.positive.list[,2] %in% dat.gene.symbol$gene_symbol,]
    true.positive.list <- t(apply(true.positive.list, 1, mixedsort))
    true.positive.list <- paste(true.positive.list[,1], true.positive.list[,2], sep = "_")

}else{
   ## Read genesets
    g.sets <- read.delim(pathways.fn, header = F, stringsAsFactors = F, row.names = 1)
    g.sets$V2 <- NULL # second column is NA
    names.gsets <- rownames(g.sets)
    g.sets <- t(g.sets)
    g.sets <- sapply(g.sets, function(x) mixedsort(strsplit(x,',')[[1]]))

    ## select genes in pathways that are also present in the set of variable genes
    genes.inpathways <- lapply(g.sets, function(x,y){
      genes.inset <- x[which(x %in% y)]
      genes.inset
      }, dat.gene.symbol$gene_symbol)
    
    unique.gene.symbols.inpathways <- unique(unlist(genes.inpathways))
    size_genes_inpathways <- length(unique.gene.symbols.inpathways)
    total_genes_symbol <- sum(!is.na(dat.gene.symbol$gene_symbol))
    print(total_genes_symbol)

    if(any(grepl("tft", pathways.fn))){
      names.gsets <- sapply(names.gsets, function(x) strsplit(x,'_')[[1]][1])
      names(genes.inpathways) <- names.gsets
      true.positive.list <- melt(genes.inpathways)
      true.positive.list <- t(apply(true.positive.list, 1, mixedsort))
      true.positive.list <- unique(paste(true.positive.list[,1], true.positive.list[,2], sep = "_"))
    }else{
      true.positive.list <- lapply(genes.inpathways, function(x){
      if(length(x) > 1){
        pathway.edge <- t(combn(mixedsort(x),2))
        pathway.edge <- paste(pathway.edge[,1], pathway.edge[,2], sep = "_")
      }else{
        pathway.edge <- NA
      }
      pathway.edge
      })
      true.positive.list <- unique(unlist(true.positive.list))
    }
    true.positive.list <- true.positive.list[which(!is.na(true.positive.list))]
}
    print("Total real edgelist size:")
    print(length(true.positive.list))

if(net_type %in% c("wgcna", "wgcna-signed", "power_wgcna")){
  # for each network within the tissue generate a fully connected subgraph
  tiss.net <- lapply(type.exp, function(x,y,z){
    print(paste("Loading network", x))
    if(net_type == "wgcna"){
    load(paste("networks", x, "wgcna_networks.Rdata", sep = "/"))
    }
    if(net_type == "wgcna-signed"){
    load(paste("networks", x, "signed_wgcna_networks.Rdata", sep = "/")) 
    }
    if(net_type == "power_wgcna"){
    load(paste("networks", x, "power_wgcna_networks.Rdata", sep = "/"))
    }
    # get module assignments for each cut-off networks
    this.net <- lapply(dat.net[[y]], function(eachNet, gene.symbol){
      eachNet <- eachNet$colors
      names(eachNet) <- gene.symbol$gene_symbol ## names for each entry in the vector -- gene names
      eachNet
    }, z)

    dat.module.network <- lapply(this.net, function(eachNet){
      total.num.modules <- mixedsort(unique(eachNet))[-1] ## remove gray modules
      if(any(table(eachNet) == 1)){
        tmp_count <- table(eachNet)
        total.num.modules <- names(tmp_count)[tmp_count != 1]
        total.num.modules <- mixedsort(as.numeric(total.num.modules))[-1]
      }
      print(total.num.modules)
      module.genes.eachNet <- lapply(total.num.modules, function(eachModule, this.network){
        genes.in.module <- names(this.network)[which(this.network == eachModule)]
        genes.in.module <- mixedsort(genes.in.module)
        if(any(is.na(genes.in.module))){
          genes.in.module <- genes.in.module[-which(is.na(genes.in.module))]
        }
        module.network <- t(combn(unique(genes.in.module), 2)) # only keep unique genes
        module.network <- paste(module.network[,1], module.network[,2], sep = "_")
        module.network
        }, eachNet)
      unlist(module.genes.eachNet)
      })
    dat.module.network
  },tiss, dat.gene.symbol)
}

if(net_type == "glasso"){
	tiss.net <- lapply(type.exp, function(x,y,z){
    print(paste("Loading dataset", y, x))
    load(paste("networks/", x, "/",y, "_glasso_networks.Rdata", sep = ""))
    

    net.edges <- lapply(dat.net, function(eachNet, z_infunc){
      if(all(sub('(^[^.]+)\\.(.*)$','\\1',colnames(eachNet),"") == z_infunc$gene_id)){
        rownames(eachNet) <- z_infunc$gene_symbol
        colnames(eachNet) <- z_infunc$gene_symbol
        na.idx <- which(is.na(rownames(eachNet)))
        eachNet <- eachNet[-na.idx,-na.idx]
        eachNet[lower.tri(eachNet)] <- t(eachNet)[lower.tri(eachNet)]
        sorted.gene.order <- order(rownames(eachNet))
        eachNet <- eachNet[sorted.gene.order, sorted.gene.order]
        eachNet[lower.tri(eachNet)] <- NA
        eachNet <- melt(eachNet, na.rm = T)
        eachNet <- paste(eachNet$Var1, eachNet$Var2, sep = "_")
        eachNet <- unique(eachNet)
    }else{
      print("Gene order and mapping error. Please manually check!")
    }
      eachNet
      }, z)
    net.edges
  },tiss, dat.gene.symbol)
}

  names(tiss.net) <- type.exp


  type.net <- vector("list", length = length(type.exp))
  names(type.net) <- type.exp 
  for(t in type.exp){
    print(t)
    type.net[[t]] <- lapply(tiss.net[[t]], function(x, z){
      total.positives <- length(x)
      total.negatives <- ((total_genes_symbol * (total_genes_symbol - 1))/2) - total.positives
      tp <- length(intersect(x,z)) # present in both network and pathways
      fn <- length(setdiff(z,x)) # present in pathways missing in network
      fp <- length(setdiff(x,z)) # present in network absent in pathways
      precision <- tp/(tp+fp)
      recall <- tp/(tp+fn)
      # c(tp, fp, fn)
      return(c(precision, recall, fn, total.positives, total.negatives))
    }, true.positive.list)
  }



    type.net.bind <- mapply(function(x,y){
      x <- do.call(rbind, x)
      colnames(x) <- c("precision", "recall", "false_negatives", "density", "total_negatives")
      x <- data.frame(x)
      x$type <- y
      x
      }, type.net, as.list(type.exp), SIMPLIFY = FALSE)

    # type.net.bind <- lapply(type.net.bind, function(x,y){
    #   rownames(x) <- y
    #   x
    #   }, threshold.values)

  # type.net.bind[[1]]$type <- "uncorrected"
  # type.net.bind[[2]]$type <- "RIN"
  # type.net.bind[[3]]$type <- "quarterPC"
  # type.net.bind[[4]]$type <- "halfPC"
  # type.net.bind[[5]]$type <- "PC"

  pr.plot <- do.call(rbind, type.net.bind)
  pr.plot$type <- factor(pr.plot$type, levels = c("pc", "halfpc", "quarterpc", "rin", "gc", "mc", "exonicRate", "raw"), labels = c("PC", "half-PC", "quarter-PC", "RIN", "gene GC%", "multi-covariate", "exonic rate", "uncorrected"))

  # plot precision and recall

  fig_pr <- ggplot(pr.plot, aes(x = recall, y = precision, col = type))+ geom_point(size = 0.3) + ggtitle(tiss) + xlim(0, 0.20)+ ylim(0.0, 0.40)
  png(paste(plot.dir, "/PR/", net_type, "_", save.fn, ".png", sep = ""), height = 2, width = 4, units = "in", res = 400)
  print(fig_pr)
  dev.off()

  ## save file
  saveRDS(pr.plot, file = paste(res.dir, "/PR/pr_density_", net_type, "_", save.fn, ".Rds", sep = "" ))
