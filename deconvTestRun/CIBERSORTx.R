
library(scater)
library(scran)



# load in processed sc data
load("sce02.RData")
# load("sce03.RData")
# load("sce04.RData")


#-----------------
 # Single cell reference file
#-----------------

sce.copy <- sce.hvg

# Relabel cells
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
    "Malignant Cells" # does not keep clustering label
  }
colData(sce.copy)$manual_annotation <- mapply(func1, sce.copy$label)


# use X02 as sc reference

# sc_ref_df <-data.frame((matrix(ncol = length(colnames(sce.copy)), nrow = length(rowData(sce.copy)$Symbol))), row.names = rowData(sce.copy)$Symbol, colnames = sce.copy$manual_annotation) # define empty data frame
# colnames(cluster_mean) <- levels(sce.hvg$label)



# convert sparse matrix to data.frame
sc_ref <- as.data.frame(as.matrix(assays(sce.copy)$counts))
# set row name to gene symbol
row.names(sc_ref) <- rowData(sce.copy)$Symbol
# set col name to cell type annotation
colnames(sc_ref) <- sce.copy$manual_annotation

# write to tab-delimited .txt
write.table(sc_ref, file = "cibersortx_sc02_ref.txt", quote = FALSE, sep = "\t", col.names = NA)




#-----------------
# Mixture file
#-----------------

library(EnsDb.Hsapiens.v86)


# load bulk data
load("gene_counts_17667X1-3.RData")  #gse

ensembl.ids <- rownames(gse)

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
ensembl.ids <- str_replace(ensembl.ids, pattern = ".[0-9]+$", replacement = "")
# convert ensembl ids to gene symbols
gene.symbols <- mapIds(EnsDb.Hsapiens.v86, keys = ensembl.ids, keytype = "GENEID", column="SYMBOL")

# build bulk reference data frame
bulk_ref <- assays(gse)$counts
row.names(bulk_ref) <- gene.symbols

# unify colnames with scRNA-seq data
colnames(bulk_ref) <- c("16030X3", "16030X2", "16030X4")

# write to tab-delimited .txt
write.table(bulk_ref, file = "cibersortx_bulk.txt", quote = FALSE, sep = "\t", col.names = NA)








