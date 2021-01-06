
# -------------------------
# scRNA-seq data processing
# -------------------------

# load sc data
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")

sce.copy <- sce

# normalize count data to library size
# Scaden author used Scanpy function “normalize_per_cell”
library(scater)
library(scuttle) # where scater functions were moved to
# Library size normalization without log transformation
sc_normalized <- normalizeCounts(sce.copy, log = FALSE) # returns a dgCMatrix (sparse)
# convert sparse matrix to data.frame
sc_normalized <- as.data.frame(as.matrix(sc_normalized))
colnames(sc_normalized) <- colnames(sce.copy)







# -----------------
# Further process in python/3.6.9
# source ~/scadenvenv/bin/activate
# -----------------

