
# Every new R session starts fresh; libraries need to be reloaded each time.

# Fix — reload the library, then rebuild pbmc by re-running Phases 1 and 2:

library(Seurat)
library(dplyr)
library(ggplot2)

# Re-run Phase 1 (data already downloaded, just reload)
pbmc.data <- Read10X(data.dir = "data/raw/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 0, min.features = 0)

# Re-run Phase 2 QC filtering
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

dim(pbmc)  # should read 13714 2638 or similar pre-doublet-filter count — check this matches what you had

# Then re-run scDblFinder:

library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, dbr = 0.023)
pbmc$scDblFinder.class <- sce$scDblFinder.class
pbmc$scDblFinder.score <- sce$scDblFinder.score

pbmc <- subset(pbmc, subset = scDblFinder.class == "singlet")
dim(pbmc)  # should read back to 2526 cells



# Normalize

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify 2,000 highly variable genes

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
top10

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


# Step 3 — Scale

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Checkpoint save (NEW)

saveRDS(pbmc, "data/processed/pbmc_phase3_normalized_scaled.rds")















