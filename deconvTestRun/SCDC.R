
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
# Note that Bisque won't run if all samples have sc and bulk, so X03 is taken out
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
pheno02.matrix <- as.data.frame(pheno02.matrix)

# use make.names() to avoid duplicate row.names
sc02.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=pheno02.matrix$Barcode, 
                       SubjectName=pheno02.matrix$Sample,
                       cellType=pheno02.matrix$manual_annotation)

sc02.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc02.eset = ExpressionSet(assayData = data.matrix(exprs02.matrix), phenoData =  new("AnnotatedDataFrame", data = sc02.pheno, varMetadata = sc02.meta))


# ------
# processing X03
# ------
sce03.copy <- sce03
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)
pheno03.matrix <- as.data.frame(pheno03.matrix)

# use make.names() to avoid duplicate row.names
sc03.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno03.matrix$Barcode, 
                         SubjectName=pheno03.matrix$Sample,
                         cellType=pheno03.matrix$manual_annotation)

sc03.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc03.eset = ExpressionSet(assayData = data.matrix(exprs03.matrix), phenoData =  new("AnnotatedDataFrame", data = sc03.pheno, varMetadata = sc03.meta))


# ------
# processing X04
# ------
sce04.copy <- sce04
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)
pheno04.matrix <- as.data.frame(pheno04.matrix)

# use make.names() to avoid duplicate row.names
sc04.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno04.matrix$Barcode, 
                         SubjectName=pheno04.matrix$Sample,
                         cellType=pheno04.matrix$manual_annotation)

sc04.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc04.eset = ExpressionSet(assayData = data.matrix(exprs04.matrix), phenoData =  new("AnnotatedDataFrame", data = sc04.pheno, varMetadata = sc04.meta))



# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs03.matrix)
sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)

# make sure the colnames are unique
colnames(sc.exprs.matrix) <- make.names(colnames(sc.exprs.matrix))

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

# use make.names() to avoid duplicate row.names
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=make.names(sc.pheno.matrix$Barcode, unique = TRUE), 
                       SubjectName=sc.pheno.matrix$Sample,
                       cellType=sc.pheno.matrix$manual_annotation)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta))



#----------------------------------------------------------------------
# RUN SCDC
# https://meichendong.github.io/SCDC/articles/SCDC.html#data-input-1
# see sample anaylsis
#----------------------------------------------------------------------
library(SCDC)
library(xbioc)

# -----
# SCDC Pre-process of scRNA-seq Data
# -----
  sc02.qc <- SCDC_qc_ONE(sc02.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X02",
                      ct.sub = c("Macrophages", "Monocytes", "Fibroblasts",  "Endothelial cells", "T-cells", "B-cells", "Malignant cells"), qcthreshold = 0.6) 

sc03.qc <- SCDC_qc_ONE(sc03.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X03",
                       ct.sub = c("Macrophages", "Monocytes", "Malignant cells"), qcthreshold = 0.6)

sc04.qc <- SCDC_qc_ONE(sc04.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X04",
                       ct.sub = c("Monocytes", "Fibroblasts",  "Endothelial cells", "Malignant cells"), qcthreshold = 0.6) 



bulk.ens <- SCDC_ENSEMBLE(bulk.eset = bulk.eset, sc.eset.list = list(X02 = sc02.qc$sc.eset.qc,  X04 = sc04.qc$sc.eset.qc), ct.varname = "cellType",
                                     sample = "SubjectName", truep = NULL, ct.sub =  c("Macrophages", "Malignant cells"), search.length = 0.01, grid.search = T)  
