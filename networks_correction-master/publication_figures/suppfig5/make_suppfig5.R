library(parallel)
library(recount)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(cowplot)
#### test association of expression PCs with gc coefficients

  
load("/work-zfs/abattle4/parsana/networks_correction/data/raw_subset.Rdata")
load("/work-zfs/abattle4/parsana/networks_correction/data/pc_loadings.Rdata")

  gc.coeff <- lapply(dat.expr, function(x) x$gc)
  # fit linear model
  gc.pc <- mapply(function(p,q,r){
    lm.gc.fit <- lm(p[,1:r] ~ q)
    lm.gc.fit
    }, pc.loadings, gc.coeff, num.pc.estimates, SIMPLIFY = FALSE
    )
  # extract p-values
  gc.pc.pvals <- lapply(gc.pc, function(lm.gc.fit){
    pval <- sapply(summary(lm.gc.fit), function(x) x$coefficients[2,4])
    pval <- p.adjust(pval, method = "BH")
    pval
    })
  # extract percent variance explained/ R2
  r2.pcs <- lapply(gc.pc, function(lm.gc.fit){
    pval <- sapply(summary(lm.gc.fit), function(x) x$r.squared)
    pval
    })

  # merge the lists to a matrix
  r2.pcs <- sapply(r2.pcs, function(x) {length(x) <- 37; x})
  gc.pc.pvals <- sapply(gc.pc.pvals, function(x) {length(x) <- 37; x})
  rownames(r2.pcs) <- paste("PC",1:37,sep="")
  rownames(gc.pc.pvals) <- paste("PC",1:37,sep="")

# prepare data for plotting with ggplot2
plot.r2 <- melt(r2.pcs)
plot.pvals <- melt(gc.pc.pvals)

  # set color pallete
  # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  myPalette <- colorRampPalette(brewer.pal(9,'Blues'), space = "Lab")

  # set default font size
  theme_set(theme_classic(base_size = 9)) 
  
  # plot p-values
  pval.fig <- ggplot(plot.pvals, aes(Var2, Var1)) + geom_tile(aes(fill = value),
       colour = "white") + 
    geom_text(size = 2, aes(label = prettyNum(value, digits=3, width=4, format="fg")
)) +
    scale_fill_gradientn(colours = rev(myPalette(10)[1:6]), limits = c(0,1)
, values = c(0, 0.005, 0.05, 0.1, 0.5, 1)
) + xlab("") + ylab("Principal Components") + theme(axis.text.x=element_text(colour="black", angle = 90), axis.text.y=element_text(colour="black"))
  ggsave("pvals.png")

  # plot r2
  r2.fig <- ggplot(plot.r2, aes(Var2, Var1)) + geom_tile(aes(fill = value),
       colour = "white") + 
    geom_text(size = 2, aes(label = prettyNum(value, digits=3, width=4, format="fg")
)) +
    scale_fill_gradientn(colours = myPalette(10)[1:6], limits = c(0,1)
    , values = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 1)
) + xlab("") + ylab("Principal Components") + theme(axis.text.x=element_text(colour="black", angle = 90), axis.text.y=element_text(colour="black"))
  ggsave("r2.png")

fig.out <- plot_grid(pval.fig + theme(legend.position="bottom", legend.direction = "horizontal"),
  r2.fig + theme(legend.position="bottom", legend.direction="horizontal"),
  align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )

pdf("suppfig5.pdf", height = 8.5, width = 7)
print(fig.out)
dev.off()

