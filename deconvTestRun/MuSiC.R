#----------------------------------------------------------------------
# Data Preparations
#----------------------------------------------------------------------


# ------
# load bulk data 
# ------
load("gene_counts_17667X1-3.RData")

# bulk sample name | single-cell sample name
# ------------------------------------------
# 17667X1          | 16030X3
# 17667X2          | 16030X2
# 17667X3          | 16030X4

# unify colnames with scRNA-seq data
colnames(gse) <- c("16030X3", "16030X2", "16030X4")

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
rownames(gse) <- str_replace(rownames(gse), pattern = ".[0-9]+$", replacement = "")

# ------
# convert bulk data to ExpressionSet 
# https://github.com/xuranw/MuSiC/issues/2
# ------
metadata <- data.frame(labelDescription= c("names"), row.names=c("names"))
# gene_exprs.matrix = assays(gse)$counts, pheno.matrix = colData(gse)
# convert DFrame to data.frame
pheno.matrix <- as.data.frame(colData(gse))
# construct ExpressionSet of bulk data
bulk.eset = ExpressionSet(assayData = data.matrix(assays(gse)$counts), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata))


# ------
# load single-cell data 
# ------
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")

# ------
# processing X02
# ------
sce.copy <- sce
# use Ensemble ID as rownames
rownames(sce.copy) <- rowData(sce.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce.copy)$counts
pheno02.matrix <- colData(sce.copy)

# ------
# processing X03
# ------
sce03.copy <- sce03
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)

# ------
# processing X04
# ------
sce04.copy <- sce04
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)


# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs03.matrix)
sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)


# ------
# combine phenotype matrix (colData)
# ------
sc.pheno.matrix <- rbind(pheno02.matrix, pheno03.matrix)
sc.pheno.matrix <- rbind(sc.pheno.matrix, pheno04.matrix)
# convert DFrame to data.frame
sc.pheno.matrix <- as.data.frame(sc.pheno.matrix)

# ------
# Generate expressionset of sc data
# https://github.com/xuranw/MuSiC/issues/2
# ------
sc.metadata <- data.frame(labelDescription= c("Sample", "Barcode", "sum", "detected", "percent_top_50", "percent_top_100", "percent_top_200", "percent_top_500", 
                                              "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent", "total", "sizeFactor", "label", "manual_annotation"), 
                       row.names=c("Sample", "Barcode", "sum", "detected", "percent_top_50", "percent_top_100", "percent_top_200", "percent_top_500", 
                                   "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent", "total", "sizeFactor", "label", "manual_annotation"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno.matrix, varMetadata = sc.metadata))




#----------------------------------------------------------------------
# RUN MuSiC
# https://xuranw.github.io/MuSiC/articles/MuSiC.html
# see sample anaylsis
#----------------------------------------------------------------------
library(MuSiC)
library(xbioc) # for pVar()
# run!
Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'manual_annotation', samples = 'Sample')



