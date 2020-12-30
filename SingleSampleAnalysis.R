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
library(tidyverse)
#Sys.setenv("DISPLAY"=":0.0")

sce <- read10xCounts(file.path("16030X2"), col.names = TRUE, type = "sparse", version = "3")


rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")
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
markers <- findMarkers(sce.hvg, groups=sce.hvg$label, pval.type="some", direction="up")


# Collect top 10 up-regulated genes from each cluster 
markers.chosen <- vector()
a <- 1:length(markers)
for (n in a) {
  marker.clust <- markers[[n]]
  markers.chosen <- append(markers.chosen, rownames(marker.clust)[1:10])
}
markers.chosen <- append(markers.chosen, "WFDC2")

plotHeatmap(sce.hvg, features=unique(markers.chosen), exprs_values = "logcounts", order_columns_by="label")



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

# subset the dataframe with selected top genes
cluster_mean_top <- subset(cluster_mean, rownames(cluster_mean) %in% unique(markers.chosen))

# format cluster_mean_top data frame from wide format to long format
dt <- cluster_mean_top %>%
       rownames_to_column() %>%
       gather(colname, value, -rowname)

# try plot the subsetted average log counts in a heatmap using ggplot2
ggplot(dt, aes(x = colname, y = rowname, fill = value)) +
  geom_tile(aes(fill = value), 
            colour = "white") +
  scale_fill_distiller(palette = "Reds", limits = c(0,10), na.value = "#de2d26",
                       direction = 1, labels = c(0.0, 2.5, 5.0, 7.5, "> 10.0"))

# use pheatmap
hm <- pheatmap(cluster_mean_top, main = "pheatmap default")

# save heatmap
save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(hm, "avgexp_heatmap1.png")

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


# ------------
# Cell Type annotation with gene expression heatmap (manual) 
# For Sample 1
# ------------
func1 <- function(x)
  if (x == "6" || x == "16") {
    "B-cells"
  } else if (x == "4" || x == "17") {
    "Macrophages"
  } else if (x == "9") {
    "Monocytes"
  } else if (x == "8" || x == "3") {
    "Fibroblasts"
  } else if (x == "7") {
    "Endothelial cells"
  } else if (x == "13") {
    "T-cells"
  } else {
    paste("Malignant cells", x) # keep clustering label
  }
colData(sce.hvg)$manual_annotation <- mapply(func1, sce.hvg$label)


# ------------
# InferCNV input preparation, see inferCNV wgene_ordering_file.txt using gene_list.txt and hg19
# Rearrage the tab-delim with command-line:
# awk '{print $4, $1, $2, $3 > "gene_ordering_file_formated.txt"}' gene_ordering_file.txt
# Get unique genes:
# awk -F , '{ a[$1]++ } END { for (b in a) { print b} }
# ------------
# Gene ordering file
write.table(data.frame(rowData(sce.hvg)$Symbol), file = "gene_list.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# get gene coordinates from UCSD hgTable as nt b > "gene_ordering_file_unique.txt"} }' gene_ordering_file.txt
df <- read.table("gene_ordering_file_formated.txt")
# eliminate duplicate genes
df <- df[!duplicated(tolower(df[,1])),]
write.table(df, file = "gene_ordering_file_unique.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



# Sample annotation file
sample_annotation <-data.frame((matrix(ncol = 2, nrow = length(colnames(sce.hvg))))) # define empty df 
sample_annotation[, 1] <- colnames(sce.hvg)
sample_annotation[, 2] <- colData(sce.hvg)$manual_annotation
write.table(sample_annotation, file = "cellAnnotations.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ------------
# Run InferCNV 
# ------------
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts(sce.hvg),
                                     annotations_file="cellAnnotations.txt",
                                     delim="\t",
                                     gene_order_file="gene_ordering_file_unique.txt",
                                     ref_group_names=c("T-cells", "B-cells", "Fibroblasts", "Monocytes", "Macrophages", "Endothelial cells"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir_keepclusterlabel",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)
