
# -------------------------
# scRNA-seq data processing
# -------------------------

# load sc data
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")

sce.copy <- sce

#---------Counts file--------------
# use X02, Include all genes


# normalize count data to library size
# Scaden author used Scanpy function “normalize_per_cell”
library(scater)
library(scuttle) # where scater functions were moved to
# Library size normalization without log transformation
sc_normalized <- normalizeCounts(sce.copy, log = FALSE) # returns a dgCMatrix (sparse)
# convert sparse matrix to data.frame
sc_normalized <- as.data.frame(as.matrix(sc_normalized))

# remove col and row names before transpose
colnames(sc_normalized) <- NULL
rownames(sc_normalized) <- NULL

# The count data should be of size (n x g), where g is the number of genes and n is the number of samples
# transpose the sc_normalized dataframe
sc_normalized <- t(sc_normalized)

# add back colnames (gene names in Ensembl ID)
colnames(sc_normalized) <- rowData(sce.copy)$ID
# add back rownames (cells, colnames of sce.copy)
rownames(sc_normalized) <- colnames(sce.copy)

# write to tab-delimited .txt
write.table(sc_normalized, file = "scaden/sc02_counts.txt", quote = FALSE, sep = "\t", col.names = NA)


# ----------------Celltype file---------------- 
# The file for the cell type labels should be of size (n x 1), 
# where n is the number of cells you have in your data. The single column in this file has to be labeled 'Celltype'.
sc_Celltype <- as.data.frame(colData(sce.copy)$manual_annotation)
colnames(sc_Celltype) <- c("Celltype")
write.table(sc_Celltype, file = "scaden/sc02_celltypes.txt", quote = FALSE, sep = "\t", row.names = FALSE)



# -------------------------
# bulk data formatting (prediction file)
# -------------------------

# load file
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
ensembl.ids <- rownames(gse)
ensembl.ids <- str_replace(ensembl.ids, pattern = ".[0-9]+$", replacement = "")
rownames(gse) <- ensembl.ids

# Generate a file of shape m X n, where m are features (genes) and n samples
# Each row corresponds to a gene, and each column to a sample
bulk <- assays(gse)$counts
write.table(bulk, file = "scaden/X020304_bulk_data.txt", quote = FALSE, sep = "\t", col.names = NA)

# -----------------
# Further deconvolution process in python/3.6.9 environment
# source ~/scadenvenv/bin/activate
# -----------------
#---------------Script-----------------
# generate training data (artificial bulk) for Scaden
# eg.  generate 1000 artificial bulk samples from 200 cells per samples with the following command
scaden simulate --cells 100 --n_samples 1000 --data scaden  --pattern *_counts.txt


# pre-process training data
scaden process data.h5ad scaden/X020304_bulk_data.txt

# model traning
scaden train processed.h5ad --model_dir scaden/model

# perform the prediction
scaden predict scaden/X020304_bulk_data.txt --model_dir scaden/model

#---------------RUN-----------------
qsub -l mem_free=20G -cwd -m e -M ychen338@jhu.edu scaden.sh