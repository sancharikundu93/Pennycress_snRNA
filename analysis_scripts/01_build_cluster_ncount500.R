#!/usr/bin/env Rscript
# Step 07.1 -- rebuild the integrated object at the nCount_RNA > 500 gate.
#
# WHY THIS EXISTS (2026-08-13): user's decision to move back to the 500 floor.
# step_06/01_build_cluster.R used nCount_RNA > 1000, which kept 18,475 nuclei and
# produced 20 clusters at res 0.6. The 500 gate keeps 25,206 (6,731 more) and gave
# 23 clusters at res 0.6 in the earlier run (jobs 652694-652720). Every current
# figure and DEG table descends from the 1000-gate object, so this rebuild starts a
# parallel lineage that has to be carried all the way through annotation before
# anything downstream can be regenerated.
#
# WHAT IS DIFFERENT FROM step_06/01 -- ONLY the nCount floor and the output path.
# Same matrices (matrices_locked), same nFeature bounds, same 2000 HVGs, same CCA
# by sample, same 30 dims, same resolutions, same seed. Keep it that way: if any
# other knob moves, the 500-vs-1000 comparison stops being attributable to the gate.
#
# WHY A NEW OUTPUT DIRECTORY. The previous 500-gate run wrote to the SAME path as
# the 1000-gate run, so building at 1000 silently destroyed the 500 object and its
# UMAPs -- ncount500_archive/ holds CSVs only, no .rds and no PDFs, which is why the
# two embeddings could not be compared without this re-run. Separate directory now,
# and the resolution UMAPs are also written as a single grid figure.
#
# Outputs the object + resolution sweep + per-resolution UMAPs + the grid.
# Annotation happens in 02 (manual), proportions in 03.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})
options(bitmapType = "cairo")
set.seed(42)

FAM <- "Times New Roman"

MAT <- "/projects/intro2gds/A_Sanchari/pennycress_snRNA_mmn30/matrices_locked"

# The gate under test. Override to reproduce the 1000-gate object for comparison:
#   NCOUNT_MIN=1000 sbatch --export=ALL run_01_build_cluster_ncount500.sh
NCOUNT_MIN <- as.numeric(Sys.getenv("NCOUNT_MIN", "500"))
NCOUNT_MAX <- 30000
NFEAT_MIN  <- 200
NFEAT_MAX  <- 8000

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
OUT  <- file.path(BASE, "ncount500_08132026/results")
FIG  <- file.path(OUT, "figures")
RDS  <- file.path(OUT, "penny_mn106_ncount500_integrated.rds")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

RESOLUTIONS <- c(0.25, 0.4, 0.6, 0.8, 1.0)

cat("gate: nFeature", NFEAT_MIN, "-", NFEAT_MAX,
    " nCount", NCOUNT_MIN, "-", NCOUNT_MAX, "\n")
cat("out :", OUT, "\n\n")

cat("########## 1. read matrices ##########\n")
samples   <- c("C1", "C2", "S1", "S2")
condition <- c(C1 = "Control", C2 = "Control", S1 = "Salt", S2 = "Salt")

objs <- lapply(samples, function(s) {
  # ReadMtx, not Read10X -- STARsolo output is uncompressed and Read10X insists on
  # gzipped barcodes.tsv.gz once it sees a v3-style features.tsv. See step_06/01.
  d <- file.path(MAT, s)
  m <- ReadMtx(mtx      = file.path(d, "matrix.mtx"),
               cells    = file.path(d, "barcodes.tsv"),
               features = file.path(d, "features.tsv"),
               feature.column = 1)
  o <- CreateSeuratObject(counts = m, project = paste0("penny_", tolower(condition[[s]])),
                          min.cells = 0, min.features = 0)
  o$sample    <- s
  o$condition <- condition[[s]]
  cat(sprintf("  %-3s %6d cells  %6d genes\n", s, ncol(o), nrow(o)))
  o
})
names(objs) <- samples

combined <- merge(objs[[1]], y = objs[-1], add.cell.ids = samples, project = "Pennycress_MN106")
cat("\nmerged:", ncol(combined), "cells\n")
print(table(combined$sample, combined$condition))

cat("\n########## 2. QC (identical gate for all four samples) ##########\n")
pre <- table(combined$sample)
penny <- subset(combined, subset = nFeature_RNA > NFEAT_MIN & nFeature_RNA < NFEAT_MAX &
                                   nCount_RNA   > NCOUNT_MIN & nCount_RNA   < NCOUNT_MAX)
post <- table(penny$sample)
qc <- data.frame(sample = names(pre), before = as.integer(pre),
                 after = as.integer(post[names(pre)]))
qc$pct_kept <- round(100 * qc$after / qc$before, 1)
qc$lost     <- qc$before - qc$after
print(qc)
cat("total retained:", ncol(penny), "\n")
write.csv(qc, file.path(OUT, "qc_cells_per_sample.csv"), row.names = FALSE)

# Per-sample UMI floor actually present after the gate. At NCOUNT_MIN=500 this is
# expected to read ~537/508/1685/691 -- i.e. the gate is nearly inert because cell
# calling already set those minima, and S1 still holds no nucleus below 1,685 while
# the other three do. That asymmetry is the known confound on any proportion claim.
floors <- tapply(penny$nCount_RNA, penny$sample, min)
meds   <- tapply(penny$nCount_RNA, penny$sample, median)
cat("\nper-sample UMI floor after gate:\n"); print(floors)
cat("per-sample median UMI after gate:\n"); print(meds)
write.csv(data.frame(sample = names(floors), floor = as.integer(floors),
                     median_umi = as.numeric(meds)),
          file.path(OUT, "per_sample_umi_floor.csv"), row.names = FALSE)

p <- VlnPlot(penny, features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample",
             ncol = 2, pt.size = 0)
ggsave(file.path(OUT, "qc_violin_post_filter.pdf"), p, width = 10, height = 5, device = cairo_pdf)

cat("\n########## 3. normalize / HVG / PCA ##########\n")
penny <- NormalizeData(penny, normalization.method = "LogNormalize", scale.factor = 10000)
penny <- FindVariableFeatures(penny, selection.method = "vst", nfeatures = 2000)
penny <- ScaleData(penny, features = rownames(penny))
penny <- RunPCA(penny, features = VariableFeatures(penny), verbose = FALSE)
ggsave(file.path(OUT, "elbow.pdf"), ElbowPlot(penny, ndims = 50),
       width = 6, height = 4, device = cairo_pdf)

cat("\n########## 4. CCA integration by sample ##########\n")
# JoinLayers FIRST -- Seurat v5 merge() leaves the object already split, and calling
# split() on it errors with "The following layers are already split".
penny[["RNA"]] <- JoinLayers(penny[["RNA"]])
penny[["RNA"]] <- split(penny[["RNA"]], f = penny$sample)
penny <- IntegrateLayers(object = penny, method = CCAIntegration,
                         orig.reduction = "pca", new.reduction = "integrated.cca",
                         verbose = FALSE)
penny[["RNA"]] <- JoinLayers(penny[["RNA"]])

penny <- FindNeighbors(penny, reduction = "integrated.cca", dims = 1:30)
penny <- RunUMAP(penny, dims = 1:30, reduction = "integrated.cca",
                 reduction.name = "umap.integrated", reduction.key = "UMAPintegrated_")

cat("\n########## 5. resolution sweep ##########\n")
rep <- data.frame()
for (r in RESOLUTIONS) {
  penny <- FindClusters(penny, resolution = r, verbose = FALSE)
  cl <- penny[[paste0("RNA_snn_res.", r)]][, 1]
  tb <- sort(table(cl), decreasing = TRUE)
  rep <- rbind(rep, data.frame(resolution = r, n_clusters = length(tb),
                               largest_cluster = names(tb)[1],
                               largest_pct = round(100 * tb[1] / sum(tb), 1),
                               second_pct  = round(100 * tb[2] / sum(tb), 1)))
  cat(sprintf("  res %.2f -> %2d clusters, largest = %.1f%%, second = %.1f%%\n",
              r, length(tb), 100 * tb[1] / sum(tb), 100 * tb[2] / sum(tb)))
}
write.csv(rep, file.path(OUT, "resolution_sweep.csv"), row.names = FALSE)

cat("\nprevious 500-gate run for reference: 17/19/23/25/25 clusters,\n")
cat("largest 19.3/19.4/12.5/12.4/10.4%. 1000-gate: 13/19/20/20/24, 26.9/16.8/14.8/13.0/10.8%.\n")

cat("\n########## 6. UMAPs ##########\n")
# Written individually AND as a grid. The grid is what survives -- the previous run's
# per-resolution PDFs were lost when the object was overwritten.
update_geom_defaults("text",  list(family = FAM))
update_geom_defaults("label", list(family = FAM))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  update_geom_defaults(ggrepel::GeomTextRepel,  list(family = FAM))
  update_geom_defaults(ggrepel::GeomLabelRepel, list(family = FAM))
}

panels <- lapply(RESOLUTIONS, function(r) {
  Idents(penny) <- paste0("RNA_snn_res.", r)
  g <- DimPlot(penny, reduction = "umap.integrated", label = TRUE, pt.size = 0.2) +
       ggtitle(paste0("resolution ", r)) + NoLegend() +
       theme(text = element_text(family = FAM),
             plot.title = element_text(family = FAM, face = "bold"))
  ggsave(file.path(OUT, sprintf("umap_res%.2f.pdf", r)), g,
         width = 7, height = 6, device = cairo_pdf)
  g
})
grid <- wrap_plots(panels, nrow = 1)
ggsave(file.path(FIG, "umap_by_resolution_ncount500.pdf"), grid,
       width = 5 * 5.2, height = 5.4, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(FIG, "umap_by_resolution_ncount500.png"), grid,
       width = 5 * 5.2, height = 5.4, dpi = 200, bg = "white", limitsize = FALSE)

for (v in c("sample", "condition")) {
  g <- DimPlot(penny, reduction = "umap.integrated", group.by = v, pt.size = 0.2) +
       theme(text = element_text(family = FAM))
  ggsave(file.path(OUT, paste0("umap_by_", v, ".pdf")), g,
         width = 7, height = 6, device = cairo_pdf)
}

saveRDS(penny, RDS)
cat("\nsaved:", RDS, "\n")
cat("\nNEXT: 02_markers_annotate.R (pick the resolution from resolution_sweep.csv first).\n")
cat("Annotation is MANUAL -- labels from the 1000-gate run do NOT carry over.\n")
cat("DONE\n")
