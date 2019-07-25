library(biomaRt)

mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                         dataset="btaurus_gene_ensembl",
                         host="http://www.ensembl.org")
go_database <- getBM(attributes = c("ensembl_gene_id",
                                    "go_id",
                                    "name_1006",
                                    "description"),
                     mart = mart )


sessionInfo()


corOptions = list(maxPOutliers =0.1)

require(matrixStats)

a = which((colMads(datExpr14_ht)) == 0)
b = which((colMads(datExpr14_cl)) == 0)
length(a);length(b)
length(unique(c(a,b)))

mad.x = apply(x,2,mad) (substitute your data for x) or colMads(x) from package matrixStats (faster than apply).
Then check which mad.x are zero.