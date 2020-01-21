load(url("http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_liver.Rdata"))

head(rse_gene)
?scale_counts()


rse <- scale_counts(rse_gene)

head((log2(rse@assays$data$counts+1)[,5]),200)

head(rse@assays$data$counts[,4])
rse_raw <- log2(rse@assays$data$counts+2)
head(rse_raw[,4])

genes_var <- apply(rse_raw, 1, var)

select_genes <- names(genes_var)[order(genes_var, decreasing = T)][1:1000]

rse_raw <- rse_raw[select_genes,]

save(rse_raw, file = "liver_subset_gtex.Rdata")



log2(2)
