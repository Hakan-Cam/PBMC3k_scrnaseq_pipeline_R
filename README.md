# PBMC3k scRNA-seq Pipeline — R/Seurat

R/Seurat implementation of the classic pbmc3k single-cell RNA-seq analysis pipeline,
built as a companion to a parallel Python/Scanpy implementation
(see: pbmc3k-scrna-pipeline).

## Pipeline
1. Environment setup (renv, Seurat 5.5.1)
2. Data loading (10x Genomics pbmc3k triplet, Read10X + CreateSeuratObject)
3. QC + doublet removal (scDblFinder) — 2,526 clean cells
4. Normalization, 2,000 HVGs, scaling
5. PCA -> UMAP -> Louvain clustering — 9 clusters
6. Marker detection (FindAllMarkers) -> 9 annotated cell types
   (CD4+ T cells, CD8+ T cells, B cells, NK cells, CD14+/FCGR3A+ Monocytes,
   Dendritic cells, Platelets)

## Notes
- Louvain used instead of Leiden due to reticulate/leidenalg environment setup issues on Windows.
- scDblFinder used instead of DoubletFinder due to Seurat v5 incompatibility.

See `progress.md` for full phase-by-phase development log.
