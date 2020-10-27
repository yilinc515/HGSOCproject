# Try analyizing 16030x2

# ------------
# Data loading
# ------------
library(ggplot)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocParallel)
library(SingleR)
library(pheatmap)
#Sys.setenv("DISPLAY"=":0.0")

sce <- read10xCounts(file.path("16030X2"), col.names = TRUE, type = "sparse", version = "3")


rownames(sce) <- uniquifyFeatureNames(
  rowData(sce)$ID, rowData(sce)$Symbol)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location


# ------------
# QC: emptyDrops() gives "no counts available to estimate the ambient profile" Error
# assume cell calling already performed, see https://support.bioconductor.org/p/123554/#123562
# ------------

unfiltered <- sce


is.mito <- grep("MT", rowData(sce)$Symbol)

stats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

sce <- addPerCellQC(sce, subsets = list(Mito=is.mito))
sce <- addPerFeatureQC(sce)
sce

#Remove worst cells
initialFilt <- sce$sum > 500 & sce$detected > 100
length(which(!initialFilt))
sce<-sce[,initialFilt]
sce



#high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
#sce <- sce[,!high.mito]



#mitochondrial read fraction and number of unique features (genes) with at least one read
metrics <- as.data.frame(colData(sce))
p <- ggplot(metrics, aes(x=detected,y=subsets_Mito_percent)) + geom_point()
p + stat_smooth(method="loess",formula=y~x,size=1,se=F,colour="blue")
p <- ggplot(metrics, aes(x=total_counts,y=subsets_Mito_percent)) + geom_point()
p <- ggplot(metrics, aes(x=detected,y=log10_total_counts)) + geom_point()

# mixture of two different distributions, the "functional" distribution and the "failed" distribution
library(flexmix)
options(scipen = 5)
set.seed(1010)
model<-flexmix(subsets_Mito_percent~detected,data=metrics,k=2)
model@components
slope_1=model@components$Comp.1[[1]]@parameters$coef[2]
intercept_1=model@components$Comp.1[[1]]@parameters$coef[1]
slope_2=model@components$Comp.2[[1]]@parameters$coef[2]
intercept_2=model@components$Comp.2[[1]]@parameters$coef[1]
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent)) + geom_point() + 
  geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) + 
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2))

# posterior probability
post <- posterior(model)
metrics$posterior_dist1 <- post[,1]
metrics$posterior_dist2 <- post[,2]
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent,colour=posterior_dist1)) + 
  geom_point() + geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) + 
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2))

# throw out cells with low posterior probability of coming from functional distribution
# cutoff = 0.25 
metrics$keep<-metrics$posterior_dist1>=0.25
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent,colour=keep)) + geom_point()
table(metrics$keep)
sce <- sce[,metrics$keep]
dim(sce)


# check data after QC
qcplots <- list()
qcplots[[1]] <- plotColData(sce, x="sum", y="subsets_Mito_percent",
                            ) + scale_x_log10()
qcplots[[2]] <-plotColData(sce, x="detected", y="subsets_Mito_percent",
                           ) + scale_x_log10()
qcplots[[3]] <-plotColData(sce, x="detected", y="sum",
                           ) + scale_y_log10()

do.call(gridExtra::grid.arrange, c(qcplots, ncol=3))
summary(high.mito)


# ------------
# Normalization
# ------------

set.seed(1000)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- logNormCounts(sce)
summary(sizeFactors(sce))

plot(librarySizeFactors(sce), sizeFactors(sce), pch=16, xlab="Library size factors", ylab="Deconvolution factors", log="xy")

# ------------
# Feature selection assuming near-Poisson variation
# ------------

# Quantifying per-gene variation
set.seed(1001)
pois <- modelGeneVarByPoisson(sce)
pois <- pois[order(pois$bio, decreasing=TRUE),]
head(pois)

plot(pois$mean, pois$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
curve(metadata(pois)$trend(x), col="dodgerblue", add=TRUE)

# Selecting highly variable genes
# select the top 10% of genes with the highest biological components
dec <- modelGeneVar(sce)
top.sce <- getTopHVGs(dec, prop=0.1)
str(top.sce)

# store hvg in AltExp
sce.hvg <- sce[top.sce,]
altExp(sce.hvg, "original") <- sce
altExpNames(sce.hvg)
# To recover original data: sce.original <- altExp(sce.hvg, "original", withColData=TRUE)


# ------------
# Dimensionality reductionn with PCA and NMF
# ------------
set.seed(1002) 
sce.hvg <- runPCA(sce.hvg) 
dim(reducedDim(sce.hvg, "PCA"))

#Choose the number of PCs Using the technical noise
set.seed(1003)
denoised.sce <- denoisePCA(sce.hvg, technical=dec)
ncol(reducedDim(denoised.sce))  # 5, elbow point gives 6

# OR: Based on population structure
# takes too long to run
#pcs <- reducedDim(sce.hvg)
#choices <- getClusteredPCs(pcs)
#metadata(choices)$chosen

# Set d = 30 for inital analysis
reducedDim(sce.hvg, "PCA") <- reducedDim(sce.hvg, "PCA")[,1:30]

# NMF  NOT WORKING!!!
set.seed(101001)
nmf.sce <- runNMF(sce.hvg, ncomponents=10, altexp = NULL)

nmf.out <- reducedDim(nmf.sce, "NMF")
nmf.basis <- attr(nmf.out, "basis")
colnames(nmf.out) <- colnames(nmf.basis) <- 1:10

per.cell <- pheatmap::pheatmap(nmf.out, silent=TRUE, 
                               main="By cell", show_rownames=FALSE,
                               color=rev(viridis::magma(100)), cluster_cols=FALSE) 

per.gene <- pheatmap::pheatmap(nmf.basis, silent=TRUE, 
                               main="By gene", cluster_cols=FALSE, show_rownames=FALSE,
                               color=rev(viridis::magma(100)))

gridExtra::grid.arrange(per.cell[[4]], per.gene[[4]], ncol=2)


# Visualization with UMAP
set.seed(1100)
sce.hvg <- runUMAP(sce.hvg, dimred="PCA")
plotReducedDim(sce.hvg, dimred="UMAP") # no available annotation



# ------------
# Clustering (Graph based)
# ------------
g <- buildSNNGraph(sce.hvg, k=20, use.dimred = 'PCA') # need to further decide what k to use
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

# colLabels(sce.hvg) <- factor(clust)  THIS IS ONLY AVAILABLE IN SCE 1.9.3, needs BioC-devel
sce.hvg$label <- factor(clust)
plotReducedDim(sce.hvg, dimred = "UMAP", colour_by="label")  
# plot by TSNE gives error "Error in `rownames<-`(`*tmp*`, value = c("AAACCCAAGCCACCGT-1", "AAACCCAAGGATGGCT-1",  : 
# attempt to set 'rownames' on an object with no dimensions"


#  Assessing cluster separation
ratio <- clusterModularity(g, clust, as.ratio=TRUE)
dim(ratio)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

# Evaluating cluster stability
myClusterFUN <- function(x) {
  g <- buildSNNGraph(x, use.dimred="PCA", type="jaccard")
  igraph::cluster_louvain(g)$membership
}

originals <- myClusterFUN(sce.hvg)

set.seed(001001)
coassign <- bootstrapCluster(sce.hvg, FUN=myClusterFUN, clusters=originals)
dim(coassign)

pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE,
         color=rev(viridis::magma(100)))






# ------------
# Marker Gene selection, one-sided t-test
# ------------

# markers <- findMarkers(sce.hvg, groups=sce.hvg$label, pval.type="some", direction="up")
markers <- findMarkers(sce.hvg, groups=sce.hvg$label, pval.type="all", direction="up")
# marker.set <- markers[["1"]]
# markers.chosen <- rownames(marker.set)[marker.set$Top <= 5]
markers.chosen <- vector()
a <- 1:length(markers)
for (n in a) {
  print(n)
  marker.clust <- markers[[n]]
  markers.chosen <- append(markers.chosen, rownames(marker.clust)[marker.clust$Top <= 10])
}

plotHeatmap(sce.hvg, features=unique(markers.chosen), order_columns_by="label")



genes <- lapply(markers, function(x) {
  rownames(x)[x$Top <= 10]
})

## uniqify the set of genes returned, after coercing to a vector
genes <- unique(unlist(genes))

plotHeatmap(sce.hvg, genes,
            colour_columns_by = "clusters",
            show_colnames = FALSE,
            clustering_method = 'ward.D2',
            fontsize_row = 6)



# initialize a df to store avg exp level for each cluster
cluster_mean <-data.frame((matrix(ncol = length(levels(sce.hvg$label)), nrow = length(row.names(sce.hvg)))), row.names = row.names(sce.hvg)) # define empty obj
colnames(cluster_mean) <- levels(sce.hvg$label)

# calculate mean logcount for each cluster
for (cluster in  levels(sce.hvg$label)) {
  
  subset <- subset(sce.hvg, , sce.hvg$label == cluster)
  
  for (feature in row.names(sce.hvg)) {
    
    cluster_mean[feature, cluster] <- mean(logcounts(subset)[feature, ])
  }
}




# ------------
# Cell Type annotation with SingleR
# ------------

# using in-built references - unsatisfactory results
ref <- BlueprintEncodeData()
pred <- SingleR(test=sce.hvg, ref=ref, labels=ref$label.main)
table(pred$labels)
plotScoreHeatmap(pred)

tab <- table(Assigned=pred$pruned.labels, Cluster=sce.hvg$label)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


# SingleR runs in two modes: (1) Single-cell: the annotation is performed for each single-cell independently
# (2) Cluster: the annotation is performed on predefined clusters, where the expression of a
# cluster is the sum expression of all cells in the cluster https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf