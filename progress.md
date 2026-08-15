---
title: "progress"
output: html_document
---

# scRNA-seq Analysis — R/Seurat Track

Companion project to `rnaseq_pipeline` (bulk RNA-seq) and `pbmc3k_scrna` (Python/Scanpy).
Goal: repeat the pbmc3k pipeline in R/Seurat, phase by phase, for a direct Scanpy-vs-Seurat comparison.

---

## 2026-08-12 — Phase 0: Environment Setup

**Status:** Complete

**Done:**
- [x] RStudio project created at `C:/Users/hcam7/PBMC3k_scrnaseq_pipeline_R`
- [x] Folder structure: `data/raw`, `data/processed`, `scripts`, `figures`, `results`
- [x] `renv` initialized for project-local package isolation
- [x] Core packages installed: Seurat 5.5.1, dplyr 1.2.1, ggplot2 4.0.3, patchwork 1.3.2, Matrix 1.7-6
- [x] DoubletFinder 2.0.6 installed via `remotes::install_github()` (later replaced — see Phase 2 notes)
- [x] Git repo initialized via `usethis::use_git()`; Git for Windows installed separately after `system("git ...")` calls failed (not on PATH)
- [x] `renv::snapshot()` run to lock package versions
- [x] Initial commit

**Notes:**
- Running natively on Windows (not WSL2) — Seurat doesn't need the bioconda toolchain that `rnaseq_pipeline` uses.
- Project was accidentally created inside `rnaseq_pipeline` on the first attempt, then a stray empty `scrnaseq_pipeline_R` folder on the second — final correct location is `PBMC3k_scrnaseq_pipeline_R`.
- Script file: `scripts/Pbmc3k scRNA Analysis_01_setup.R`

---

## 2026-08-13 — Phase 1: Data Loading

**Status:** Complete

**Done:**
- [x] Downloaded pbmc3k triplet (barcodes.tsv, genes.tsv, matrix.mtx) from 10x Genomics
- [x] Loaded raw/unfiltered into Seurat object via `Read10X()` + `CreateSeuratObject()` (min.cells=0, min.features=0)
- [x] Confirmed `dim(pbmc)` = 32738 x 2700 (matches Python AnnData 2,700 x 32,738, transposed)

**Notes:**
- Seurat convention: features (genes) as rows, cells as columns — opposite orientation from AnnData, same data.
- Script file: `scripts/Pbmc3k scRNA Analysis_02_load_data.R`

---

## 2026-08-13 — Phase 2: QC + Doublet Removal

**Status:** Complete

**Done:**
- [x] Calculated percent.mt via `PercentageFeatureSet()`
- [x] Visualized QC metrics (VlnPlot: nFeature_RNA, nCount_RNA, percent.mt) — matched expected pbmc3k distribution shape
- [x] Applied standard Seurat tutorial thresholds: nFeature_RNA 200-2500, percent.mt < 5
- [x] Attempted DoubletFinder — failed with persistent `cannot xtfrm data frames` error (Seurat v5 assay incompatibility), even after `Assay` class conversion and full pipeline rebuild
- [x] Switched to scDblFinder (Bioconductor, actively maintained, native Seurat v5 support) — 2,526 singlets / 112 doublets detected (dbr=0.023)
- [x] Filtered to singlets only — final clean cell count: **2,526**

**Notes:**
- Python/Scrublet clean count: 2,595. R/scDblFinder clean count: 2,526 — ~2.6% difference, expected given different doublet-detection algorithms, not an error.
- DoubletFinder is not currently reliable with Seurat v5 objects as of this session (Aug 2026) — scDblFinder recommended as the modern alternative going forward.
- Git housekeeping issue caught and fixed this phase: `data/raw/` was never actually in `.gitignore` (only R-boilerplate ignores from `usethis` defaults were present), so raw data files got committed. Fixed via `git rm -r --cached data/raw` + adding `data/raw/` and `renv/library/` to `.gitignore`.

---

## 2026-08-13 — Phase 3: Normalization + HVGs + Scaling

**Status:** In progress

**Done:**
- [ ] Log-normalization via `NormalizeData()`
- [ ] 2,000 HVGs identified via `FindVariableFeatures(selection.method="vst")`
- [ ] Data scaled via `ScaleData()`

**Notes:**
- A preliminary normalize/HVG/scale/PCA pass was already run once in Phase 2 solely to feed scDblFinder's PCA-space requirement — this phase redoes it properly on the clean, filtered 2,526-cell object for the record.

**Next session:** Continue Phase 3 (awaiting VariableFeaturePlot + ScaleData confirmation), then Phase 4 — PCA → UMAP → Leiden/Louvain clustering.

## 2026-08-14 — Phase 3: Normalization + HVGs + Scaling

**Status:** Complete

**Done:**
- [x] Log-normalization via NormalizeData() (LogNormalize, scale.factor=10000)
- [x] 2,000 HVGs identified via FindVariableFeatures(selection.method="vst")
- [x] Top 10 HVGs confirmed biologically sensible: PPBP, LYZ, S100A9, FTL, IGLL5, GNLY, FTH1, PF4, HLA-DRA, S100A8
- [x] Data scaled via ScaleData() across all genes
- [x] Checkpoint saved: data/processed/pbmc_phase3_normalized_scaled.rds

**Notes:**
- Started new R session this phase with empty Environment — pbmc object had to be fully rebuilt from Phase 1+2 steps since nothing was saved to disk previously. Introduced .rds checkpoint saving from this phase forward to prevent repeat rework.
- data/processed/*.rds added to .gitignore — these are reproducible from scripts, no need to version large binary files in Git.
- VariableFeaturePlot() threw a harmless "log-10 transformation introduced infinite values" warning — expected due to near-zero-expression genes on the log10 x-axis, does not affect results.

**Next session:** Phase 4 — PCA → UMAP → clustering (Louvain default, or Leiden via algorithm=4 for closer match to Python).
## 2026-08-14 — Phase 4: PCA + Clustering + UMAP

**Status:** Complete

**Done:**
- [x] PCA via RunPCA() on 2,000 HVGs, checked ElbowPlot
- [x] FindNeighbors(dims=1:10) + FindClusters(resolution=0.5) — Louvain algorithm (Seurat default)
- [x] 9 clusters identified (0-8) — same cluster count as Python/Leiden run
- [x] UMAP via RunUMAP(dims=1:10) — visually consistent structure with Python output (large connected lymphocyte mass, isolated B/NK-like cluster, separate monocyte-like group)
- [x] Checkpoint saved: data/processed/pbmc_phase4_clustered.rds

**Notes:**
- Attempted Leiden clustering (algorithm=4) to exactly match Python, but reticulate/leidenalg environment setup (conda env, numpy/leidenalg install) hit repeated snags and was abandoned in favor of Seurat's default Louvain algorithm.
- Despite different clustering algorithms (Louvain vs Python's Leiden), cluster count (9) and UMAP topology are consistent — good sign of underlying biological signal robustness across methods.

**Next session:** Phase 5 — Marker detection via FindAllMarkers(), biological annotation of all 9 clusters.
## 2026-08-15 — Phase 5: Marker Detection + Cell Type Annotation

**Status:** Complete

**Done:**
- [x] FindAllMarkers() run across all 9 clusters (Wilcoxon, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
- [x] Top 5 markers per cluster reviewed
- [x] DotPlot + FeaturePlot of canonical markers (CD3D, MS4A1, CD14, GNLY, PPBP, FCGR3A, etc.) confirmed cluster identities
- [x] All 9 clusters named: CD4+ T cells (memory), CD14+ Monocytes, CD4+ T cells (naive), B cells, CD8+ T cells, FCGR3A+ Monocytes, NK cells, Dendritic cells, Platelets
- [x] RenameIdents() applied, labeled UMAP generated
- [x] Final checkpoint saved: data/processed/pbmc_phase5_annotated.rds

**Notes:**
- All 9 canonical pbmc3k cell types recovered, matching Python/Scanpy pipeline's cell type count exactly.
- Marker-based identity calls were unambiguous — each cluster had a clear, high-confidence canonical marker signature.

---

## R/Seurat Pipeline — COMPLETE

Full pipeline recap:
- [x] Phase 0 — Environment (PBMC3k_scrnaseq_pipeline_R), renv, Git
- [x] Phase 1 — pbmc3k triplet downloaded, loaded into Seurat (32,738 x 2,700)
- [x] Phase 2 — QC + scDblFinder doublet removal -> 2,526 clean cells (Python: 2,595)
- [x] Phase 3 — Normalization, 2,000 HVGs, scaling
- [x] Phase 4 — PCA -> UMAP -> Louvain (9 clusters; Leiden attempted, abandoned due to reticulate setup issues)
- [x] Phase 5 — Marker detection -> 9 named, biologically validated cell types
- [x] Committed to Git throughout

**Python vs R comparison:**
| | Python/Scanpy | R/Seurat |
|---|---|---|
| Clean cells | 2,595 (Scrublet) | 2,526 (scDblFinder) |
| Clustering algorithm | Leiden | Louvain |
| Clusters found | 9 | 9 |
| Cell types identified | 9 | 9 |






































