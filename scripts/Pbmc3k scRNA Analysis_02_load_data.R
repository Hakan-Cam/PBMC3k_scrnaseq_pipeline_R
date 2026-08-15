getwd()


# Download the data
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
destfile <- "data/raw/pbmc3k_filtered_gene_bc_matrices.tar.gz"

download.file(url, destfile, mode = "wb")
untar(destfile, exdir = "data/raw")

# Verify

list.files("data/raw/filtered_gene_bc_matrices/hg19/")

# Load into Seurat, unfiltered

library(Seurat)

pbmc.data <- Read10X(data.dir = "data/raw/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 0,
  min.features = 0
)

pbmc
dim(pbmc)

# Calculate %MT and inspect QC metrics

library(Seurat)
library(dplyr)
library(ggplot2)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize distributions before deciding on cutoffs

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Apply standard Seurat tutorial thresholds

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

dim(pbmc)

# Doublet detection with scDblFinder 

# Install

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scDblFinder")
BiocManager::install("SingleCellExperiment")

# Run
library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(pbmc)

sce <- scDblFinder(sce, dbr = 0.023)  # expected doublet rate for ~2,700 cells

pbmc$scDblFinder.class <- sce$scDblFinder.class
pbmc$scDblFinder.score <- sce$scDblFinder.score

table(pbmc$scDblFinder.class)

file.exists("progress.md")









