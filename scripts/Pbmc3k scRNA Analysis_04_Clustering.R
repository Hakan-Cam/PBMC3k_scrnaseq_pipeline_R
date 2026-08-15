#P hase 4 — PCA → Neighbors → Louvain clustering → UMAP

# Step 1 — Reload checkpoint

library(Seurat)
pbmc <- readRDS("data/processed/pbmc_phase3_normalized_scaled.rds")


# Step 2 — PCA

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

ElbowPlot(pbmc)

# Step 3 — Neighbors + Louvain clustering

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)  # algorithm=1 (Louvain) is the default, no need to specify

table(Idents(pbmc))

# Step 4 — UMAP

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Step 5 — Checkpoint save

saveRDS(pbmc, "data/processed/pbmc_phase4_clustered.rds")











