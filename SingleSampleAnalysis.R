# Try analyizing 16030x2

# ------------
# Data loading
# ------------
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocParallel)
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


'''
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- high.mito


qcplots <- list()
qcplots[[1]] <- plotColData(unfiltered, x="sum", y="subsets_Mito_percent",
                            colour_by="discard") + scale_x_log10()
qcplots[[2]] <-plotColData(unfiltered, x="detected", y="subsets_Mito_percent",
                         colour_by="discard") + scale_x_log10()
qcplots[[3]] <-plotColData(unfiltered, x="detected", y="sum",
                           colour_by="discard") + scale_y_log10()

do.call(gridExtra::grid.arrange, c(qcplots, ncol=3))
summary(high.mito)

metrics <- as.data.frame(colData(sce))
'''

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
plot(librarySizeFactors(sce), sizeFactors(sce), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")

# ------------
# Feature selection
# ------------

