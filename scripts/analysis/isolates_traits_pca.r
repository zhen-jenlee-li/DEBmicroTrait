library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggplot2)

names <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_name.csv") 
dat <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_pca.csv") 
rownames(dat) <- names[,1]
colnames(dat) <- c("Genome size", "Minimum generation time", "RRN", "Total transporters", "Sugars", "Auxins", "Fatty acids", "Nucleotides", "Amino acids", "Organic acids", "GHs", "Response")
pca.res <- PCA(dat, quali.sup = c(12), scale.unit = TRUE, graph = FALSE)

p <- fviz_pca_biplot(pca.res,
             geom.ind = c("text"), labelsize=3, # show points only (nbut not "text")
             col.ind = dat$Response, # color by groups
             palette = c("#bd5e8b", "#007100", '#525252'),
             addEllipses = FALSE, # Concentration ellipses
             pointshape = 15, pointsize = 3,
             geom.var = c("arrow", "point", "text"), col.var="black", alpha.var=0.3, repel=TRUE,
             title="", show.legend=FALSE, xlab="PC1 (31% contribution to total variance)", ylab="PC2 (17.7% contribution to total variance)")

p <- fviz_add(p, pca.res$quali.sup$coord, color = c("#bd5e8b", "#007100", '#525252'), labelsize=5, pointsize=1, repel=FALSE)
ggsave("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/plots/IsolatesPCA.png", width=5.5, height=5.5)
#ggsave("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/plots/IsolatesPCA.eps", width=5.5, height=5.5)
