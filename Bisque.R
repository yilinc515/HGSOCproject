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
# load("sce03.RData")
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
# X02+X04
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs04.matrix)

# make sure the colnames are unique
colnames(sc.exprs.matrix) <- make.names(colnames(sc.exprs.matrix))

# ------
# combine phenotype matrix (colData)
# ------
# X02+X04
sc.pheno.matrix <- rbind(pheno02.matrix, pheno04.matrix)

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
# Run Bisque deconvolution
# Use Reference-based decomposition mode
#----------------------------------------------------------------------
library(Biobase)
library(BisqueRNA)

# without markers: By default, Bisque uses all genes for decomposition.
# may supply a list of genes (such as marker genes) to be used with the markers parameter, see vignette
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)
# A list is returned with decomposition estimates in slot bulk.props.

ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)