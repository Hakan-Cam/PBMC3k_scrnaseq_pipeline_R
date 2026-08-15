
# Load Phase 4's detection

pbmc <- readRDS("data/processed/pbmc_phase4_clustered.rds")

# Load needed Libraries

library(Seurat)
library(dplyr)

# Find markers for every cluster

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 5 markers per cluster, ranked by fold change
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

print(top5, n = 50)

# Visualize key canonical markers

canonical_markers <- c(
  "IL7R", "CD3D",        # T cells (general)
  "CD8A",                # CD8 T cells
  "MS4A1",                # B cells
  "GNLY", "NKG7",         # NK cells
  "CD14", "LYZ",          # CD14+ Monocytes
  "FCGR3A", "MS4A7",      # FCGR3A+ Monocytes
  "FCER1A", "CST3",       # Dendritic cells
  "PPBP"                  # Platelets
)

DotPlot(pbmc, features = canonical_markers) + RotatedAxis()

# Featureplot for a quick visual gut-check (optional but helpful)

FeaturePlot(pbmc, features = c("CD3D", "MS4A1", "CD14", "GNLY", "PPBP", "FCGR3A"))

# Apply labels to the object

new.cluster.ids <- c(
  "CD4+ T cells (memory)",   # 0
  "CD14+ Monocytes",          # 1
  "CD4+ T cells (naive)",     # 2
  "B cells",                  # 3
  "CD8+ T cells",              # 4
  "FCGR3A+ Monocytes",         # 5
  "NK cells",                  # 6
  "Dendritic cells",           # 7
  "Platelets"                  # 8
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()


# Save in Figure File 

library(ggplot2)

# Recreate the plot object (or reuse it if still in your environment as a saved variable)

final_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()


ggsave(
  filename = "figures/phase5_annotated_umap.png",
  plot = final_umap,
  width = 7, height = 6, dpi = 300
)

















