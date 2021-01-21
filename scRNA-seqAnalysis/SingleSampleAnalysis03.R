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

sce03 <- read10xCounts(file.path("16030X3"), col.names = TRUE, type = "sparse", version = "3")


rownames(sce03) <- uniquifyFeatureNames(rowData(sce03)$ID, rowData(sce03)$Symbol)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce03)$ID, column="SEQNAME", keytype="GENEID")
rowData(sce03)$CHR <- location


# ------------
# QC: emptyDrops() gives "no counts available to estimate the ambient profile" Error
# assume cell calling already performed, see https://support.bioconductor.org/p/123554/#123562
# ------------

unfiltered <- sce03


is.mito <- grep("MT", rowData(sce03)$Symbol)

stats <- perCellQCMetrics(sce03, subsets=list(Mito=is.mito))

sce03 <- addPerCellQC(sce03, subsets = list(Mito=is.mito))
sce03 <- addPerFeatureQC(sce03)
sce03

#Remove worst cells
initialFilt <- sce03$sum > 500 & sce03$detected > 100
length(which(!initialFilt))
sce03<-sce03[,initialFilt]
sce03



#high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
#sce03 <- sce03[,!high.mito]



#mitochondrial read fraction and number of unique features (genes) with at least one read
metrics <- as.data.frame(colData(sce03))
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
sce03 <- sce03[,metrics$keep]
dim(sce03)


# check data after QC
qcplots <- list()
qcplots[[1]] <- plotColData(sce03, x="sum", y="subsets_Mito_percent",
                            ) + scale_x_log10()
qcplots[[2]] <-plotColData(sce03, x="detected", y="subsets_Mito_percent",
                           ) + scale_x_log10()
qcplots[[3]] <-plotColData(sce03, x="detected", y="sum",
                           ) + scale_y_log10()

do.call(gridExtra::grid.arrange, c(qcplots, ncol=3))
summary(high.mito)


# ------------
# Normalization
# ------------

set.seed(1000)
clusters <- quickCluster(sce03)
sce03 <- computeSumFactors(sce03, cluster=clusters)
sce03 <- logNormCounts(sce03)
summary(sizeFactors(sce03))

plot(librarySizeFactors(sce03), sizeFactors(sce03), pch=16, xlab="Library size factors", ylab="Deconvolution factors", log="xy")

# ------------
# Feature selection assuming near-Poisson variation
# ------------

# Quantifying per-gene variation
set.seed(1001)
pois <- modelGeneVarByPoisson(sce03)
pois <- pois[order(pois$bio, decreasing=TRUE),]
head(pois)

plot(pois$mean, pois$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
curve(metadata(pois)$trend(x), col="dodgerblue", add=TRUE)

# Selecting highly variable genes
# select the top 10% of genes with the highest biological components
dec <- modelGeneVar(sce03)
top.sce03 <- getTopHVGs(dec, prop=0.1)
str(top.sce03)

# store hvg in AltExp
sce03.hvg <- sce03[top.sce03,]
altExp(sce03.hvg, "original") <- sce03
altExpNames(sce03.hvg)
# To recover original data: sce03.original <- altExp(sce03.hvg, "original", withColData=TRUE)


# ------------
# Dimensionality reductionn with PCA and NMF
# ------------
set.seed(1002) 
sce03.hvg <- runPCA(sce03.hvg) 
dim(reducedDim(sce03.hvg, "PCA"))

#Choose the number of PCs Using the technical noise
set.seed(1003)
denoised.sce03 <- denoisePCA(sce03.hvg, technical=dec)
ncol(reducedDim(denoised.sce03))  # 5, elbow point gives 6

# OR: Based on population structure
# takes too long to run
#pcs <- reducedDim(sce03.hvg)
#choices <- getClusteredPCs(pcs)
#metadata(choices)$chosen

# Set d = 30 for inital analysis
reducedDim(sce03.hvg, "PCA") <- reducedDim(sce03.hvg, "PCA")[,1:30]

# NMF  NOT WORKING!!!
set.seed(101001)
nmf.sce03 <- runNMF(sce03.hvg, ncomponents=10, altexp = NULL)

nmf.out <- reducedDim(nmf.sce03, "NMF")
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
sce03.hvg <- runUMAP(sce03.hvg, dimred="PCA")
plotReducedDim(sce03.hvg, dimred="UMAP") # no available annotation



# ------------
# Clustering (Graph based)
# ------------
g03 <- buildSNNGraph(sce03.hvg, k=10, use.dimred = 'PCA') # need to further decide what k to use  # change k = 10 to get a finer resolution
clust03 <- igraph::cluster_walktrap(g03)$membership
table(clust03)

# colLabels(sce03.hvg) <- factor(clust)  THIS IS ONLY AVAILABLE IN sce03 1.9.3, needs BioC-devel
sce03$label <- factor(clust03)
sce03.hvg$label <- factor(clust03)

plotReducedDim(sce03.hvg, dimred = "UMAP", colour_by="label")  
# plot by TSNE gives error "Error in `rownames<-`(`*tmp*`, value = c("AAACCCAAGCCACCGT-1", "AAACCCAAGGATGGCT-1",  : 
# attempt to set 'rownames' on an object with no dimensions"


#  Assessing cluster separation
ratio <- clusterModularity(g03, clust03, as.ratio=TRUE)
dim(ratio)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

# Evaluating cluster stability
myClusterFUN <- function(x) {
  g03 <- buildSNNGraph(x, use.dimred="PCA", type="jaccard")
  igraph::cluster_louvain(g03)$membership
}

originals <- myClusterFUN(sce03.hvg)

set.seed(001001)
coassign <- bootstrapCluster(sce03.hvg, FUN=myClusterFUN, clusters=originals)
dim(coassign)

pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE,
         color=rev(viridis::magma(100)))






# ------------
# Marker Gene selection, one-sided t-test
# ------------

# markers <- findMarkers(sce03.hvg, groups=sce03.hvg$label, pval.type="some", direction="up")
markers <- findMarkers(sce03.hvg, groups=sce03.hvg$label, pval.type="some", direction="up")


# Collect top 10 up-regulated genes from each cluster 
markers.chosen <- vector()
a <- 1:length(markers)
for (n in a) {
  marker.clust <- markers[[n]]
  markers.chosen <- append(markers.chosen, rownames(marker.clust)[1:10])
}
markers.chosen <- append(markers.chosen, "WFDC2")

plotHeatmap(sce03.hvg, features=unique(markers.chosen), exprs_values = "logcounts", order_columns_by="label")



genes <- lapply(markers, function(x) {
  rownames(x)[x$Top <= 10]
})

## uniqify the set of genes returned, after coercing to a vector
genes <- unique(unlist(genes))

plotHeatmap(sce03.hvg, genes,
            colour_columns_by = "clusters",
            show_colnames = FALSE,
            clustering_method = 'ward.D2',
            fontsize_row = 6)



# initialize a df to store avg exp level for each cluster
cluster_mean <-data.frame((matrix(ncol = length(levels(sce03.hvg$label)), nrow = length(row.names(sce03.hvg)))), row.names = row.names(sce03.hvg)) # define empty obj
colnames(cluster_mean) <- levels(sce03.hvg$label)

# calculate mean logcount for each cluster
for (cluster in  levels(sce03.hvg$label)) {
  
  subset <- subset(sce03.hvg, , sce03.hvg$label == cluster)
  
  for (feature in row.names(sce03.hvg)) {
    
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
pred <- SingleR(test=sce03.hvg, ref=ref, labels=ref$label.main)
table(pred$labels)
plotScoreHeatmap(pred)

tab <- table(Assigned=pred$pruned.labels, Cluster=sce03.hvg$label)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


# ------------
# Cell Type annotation with gene expression heatmap (manual) 
# For Sample 1
# ------------
func1 <- function(x)
  if (x == "3" || x == "9") {
    "Monocytes"
  } else if (x == "8" || x == "10") {
    "Macrophages"
  } else if (x == "5" || x == "6" || x == "7") {
    "T-cells"
  } else if (x == "1") {
    "Ambiguous"
  } else {
    # paste("Malignant cells", x) # keep clustering label for inferCNV
    "Malignant Cells"
  }
colData(sce03)$manual_annotation <- mapply(func1, sce03$label)
colData(sce03.hvg)$manual_annotation <- mapply(func1, sce03.hvg$label)


#-------------
# Save data
#-------------
save(sce03.hvg, sce03, file = "sce03.RData")



# ------------
# InferCNV input preparation, see inferCNV wgene_ordering_file.txt using gene_list03.txt and hg19
# Rearrage the tab-delim with command-line:
# awk '{print $4, $1, $2, $3 > "gene_ordering_file_formated03.txt"}' gene_ordering_file03.txt
# Get unique genes:
# awk -F , '{ a[$1]++ } END { for (b in a) { print b} }
# ------------
# Gene ordering file
write.table(data.frame(rowData(sce03.hvg)$Symbol), file = "gene_list03.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# get gene coordinates from UCSD hgTable as nt b > "gene_ordering_file_unique03.txt"} }' gene_ordering_file_formated03.txt
df <- read.table("gene_ordering_file_formated03.txt")
# eliminate duplicate genes
df <- df[!duplicated(tolower(df[,1])),]
write.table(df, file = "gene_ordering_file_unique03.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



# Sample annotation file
sample_annotation <-data.frame((matrix(ncol = 2, nrow = length(colnames(sce03.hvg))))) # define empty df 
sample_annotation[, 1] <- colnames(sce03.hvg)
sample_annotation[, 2] <- colData(sce03.hvg)$manual_annotation
write.table(sample_annotation, file = "cellAnnotations03.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ------------
# Run InferCNV 
# ------------
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts(sce03.hvg),
                                     annotations_file="cellAnnotations03.txt",
                                     delim="\t",
                                     gene_order_file="gene_ordering_file_unique03.txt",
                                     ref_group_names=c("T-cells", "Monocytes", "Macrophages"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir_keepclusterlabel_03",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)

# Convert to Seurat

sce03_seurat <- Convert(from = sce03.hvg, to = "seurat")
