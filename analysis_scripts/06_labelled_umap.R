#!/usr/bin/env Rscript
# Step 06 -- the annotated UMAP for the 500-gate object.
#
# Four figures, because they answer different questions:
#   a  numbered      -- cluster IDs, so every other table can be traced back
#   b  labelled      -- cell-type names on the plot
#   c  split         -- Control | Salt side by side, same embedding
#   d  numbered+key  -- IDs on the plot with a legend giving ID -> name, which is the
#                       version that stays readable when several clusters share a
#                       label (three Endodermis, two Atrichoblast, two Pericycle,
#                       two S phase, two Cortex here -- names alone collide)
#
# Labels come from 02_annotation/labels_res0.6.csv. Provenance of that file:
#   - 19 clusters settled by >=2 of 3 independent lines (crosstab against the
#     curated 1000-gate labels, root-atlas hypergeometric, OMG)
#   - clusters 4, 5, 15 set by the user 2026-08-13
#   - clusters 13, 14 taken from the crosstab majority (74.9% / 66.2%)
#
# Outputs -> ncount500_08132026/results/06_labelled_umap/

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel)
})
options(bitmapType = "cairo")

FAM  <- "Times New Roman"
BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
LAB  <- file.path(ROOT, "results/02_annotation/labels_res0.6.csv")
OUT  <- file.path(ROOT, "results/06_labelled_umap")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CLUSTER_COL <- Sys.getenv("CLUSTER_COL", "RNA_snn_res.0.6")
REDUC <- "umap.integrated"
COND  <- c(Control = "#F0BC35", Salt = "#C2185B")

## ---- cell-type palette ------------------------------------------------------
# Single shared definition, so the UMAP and the proportion bars cannot drift apart.
source(file.path(ROOT, "scripts/palette_celltypes.R"))

stopifnot(file.exists(RDS), file.exists(LAB))

# DimPlot's on-plot labels are drawn by geom_text/ggrepel, which do NOT inherit the
# theme font -- their defaults have to be set directly or the labels come out sans.
update_geom_defaults("text",  list(family = FAM))
update_geom_defaults("label", list(family = FAM))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  update_geom_defaults(ggrepel::GeomTextRepel,  list(family = FAM))
  update_geom_defaults(ggrepel::GeomLabelRepel, list(family = FAM))
}

obj <- readRDS(RDS)
stopifnot(CLUSTER_COL %in% colnames(obj@meta.data), REDUC %in% Reductions(obj))
# shared orientation (rotate_90) -- see palette_celltypes.R
obj <- orient_reduction(obj, REDUC)
cat("orientation:", ORIENT, "\n")

lab <- fread(LAB, colClasses = "character")
stopifnot(all(c("cluster", "cell_type") %in% names(lab)))
cl  <- as.character(obj@meta.data[[CLUSTER_COL]])
hit <- match(cl, lab$cluster)
if (anyNA(hit)) {
  # never let an unlabelled cluster vanish -- that is how a cell type silently
  # disappears from every downstream table
  miss <- sort(unique(cl[is.na(hit)]))
  stop("no label for cluster(s): ", paste(miss, collapse = ", "))
}
obj$cell_type <- lab$cell_type[hit]
obj$condition <- factor(as.character(obj$condition), levels = names(COND))

## ---- resolve colours --------------------------------------------------------
# Fail loudly on an unknown cell type rather than let ggplot drop those nuclei to
# grey -- a silently grey cluster is indistinguishable from a real "unassigned".
present <- sort(unique(obj$cell_type))
unknown <- setdiff(present, names(PAL))
if (length(unknown)) stop("no colour in PAL for: ", paste(unknown, collapse = ", "))
obj$cell_type <- factor(obj$cell_type, levels = intersect(CT_LEVELS, present))
stopifnot(!anyNA(obj$cell_type))
CT_COLS <- PAL[levels(obj$cell_type)]

# per-cluster colours and _1/_2/_3 suffixing: cluster_label_map() in
# palette_celltypes.R. Moved there 2026-08-14 unchanged, so 38_umap_split_by_replicate.R
# can use the identical labels rather than restating the logic and drifting.
clm         <- cluster_label_map(cl, as.character(obj$cell_type))
map         <- clm$map
CL_COLS     <- clm$cl_cols
ctn_levels  <- clm$ctn_levels
CTN_COLS    <- clm$ctn_cols
obj$cell_type_num <- factor(clm$label_per_cell, levels = ctn_levels)
stopifnot(!anyNA(obj$cell_type_num))
cat("\ncluster -> colour:\n"); print(map[order(as.numeric(cluster)),
                                         .(cluster, cell_type, colour)])

n  <- table(obj$cell_type)
cat("nuclei:", ncol(obj), " cell types:", length(n), "\n")
print(sort(n, decreasing = TRUE))

sz <- as.data.table(obj@meta.data)[, .(cells = .N), by = .(cluster = get(CLUSTER_COL),
                                                           cell_type)]
sz[, pct := round(100 * cells / sum(cells), 2)]
setorder(sz, -cells)
fwrite(sz, file.path(OUT, "cluster_label_sizes.csv"))

# Type sizes. The previous run set base_size 12 with 9 pt legends and 2.8-4 pt on-plot
# labels, which is unreadable once the figure is scaled into a document -- the whole
# point of these panels is that a reader can name the clusters. Everything below is
# sized to stay legible at print width; the canvases grow with it so nothing crowds.
th <- theme_classic(base_family = FAM, base_size = 16) +
  theme(text = element_text(colour = "black", family = FAM),
        axis.text = element_text(colour = "black", size = 13),
        axis.title = element_text(colour = "black", size = 15),
        plot.title = element_text(family = FAM, face = "bold", size = 18),
        legend.text = element_text(family = FAM, size = 13),
        legend.title = element_text(family = FAM, size = 14),
        strip.text = element_text(family = FAM, colour = "black", size = 15))

## ---- label positions, on the ORIENTED embedding -----------------------------
# Centroids must be computed after orient_reduction(), or the labels land where the
# clusters used to be.
emb <- as.data.table(Embeddings(obj, REDUC)[, 1:2])
setnames(emb, c("x", "y"))
emb[, cluster := cl]
emb[, lbl := as.character(obj$cell_type_num)]
cen_cluster <- emb[, .(x = median(x), y = median(y)), by = .(label = cluster)]
cen_label   <- emb[, .(x = median(x), y = median(y)), by = .(label = lbl)]

## ---- a. numbered ------------------------------------------------------------
pa <- add_labels(
  DimPlot(obj, reduction = REDUC, group.by = CLUSTER_COL, label = FALSE,
          pt.size = 0.2) +
    scale_colour_manual(values = CL_COLS) +
    ggtitle(sprintf("Clusters (res 0.6, %s nuclei)", format(ncol(obj), big.mark = ","))) +
    th + NoLegend(),
  cen_cluster)
ggsave(file.path(OUT, "Fig_umap_numbered.pdf"), pa, device = cairo_pdf,
       width = 9, height = 9)
ggsave(file.path(OUT, "Fig_umap_numbered.png"), pa, width = 9, height = 9,
       dpi = 200, bg = "white")

## ---- b. labelled ------------------------------------------------------------
pb <- add_labels(
  DimPlot(obj, reduction = REDUC, group.by = "cell_type_num", label = FALSE,
          pt.size = 0.2) +
    scale_colour_manual(values = CTN_COLS) +
    ggtitle("Annotated cell types") + th + NoLegend(),
  cen_label)
ggsave(file.path(OUT, "Fig_umap_labelled.pdf"), pb, device = cairo_pdf,
       width = 9.5, height = 9.5)
ggsave(file.path(OUT, "Fig_umap_labelled.png"), pb, width = 9.5, height = 9.5,
       dpi = 200, bg = "white")

# same panel keeping the legend as well
pb2 <- add_labels(
  DimPlot(obj, reduction = REDUC, group.by = "cell_type_num", label = FALSE,
          pt.size = 0.2) +
    scale_colour_manual(values = CTN_COLS, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    ggtitle("Annotated cell types") + th,
  cen_label)
ggsave(file.path(OUT, "Fig_umap_labelled_legend.pdf"), pb2, device = cairo_pdf,
       width = 12.5, height = 9.5)
ggsave(file.path(OUT, "Fig_umap_labelled_legend.png"), pb2, width = 12.5, height = 9.5,
       dpi = 200, bg = "white")

## ---- c. Control | Salt ------------------------------------------------------
# Names BOTH on the plot and in the legend. Centroids come from ALL nuclei rather
# than per facet, so Control and Salt carry identical label positions and the eye can
# compare the same point across the two panels.
# The Control/Salt strip is left unboxed -- theme_classic draws a bordered rect around
# facet titles by default -- and set bold instead.
pc <- add_labels(
  DimPlot(obj, reduction = REDUC, group.by = "cell_type_num",
          split.by = "condition", label = FALSE, pt.size = 0.2) +
    scale_colour_manual(values = CTN_COLS, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    ggtitle("Annotated cell types, Control vs Salt") + th +
    theme(strip.background = element_blank(),
          strip.text = element_text(family = FAM, colour = "black", size = 17,
                                    face = "bold")),
  cen_label)
ggsave(file.path(OUT, "Fig_umap_labelled_by_condition.pdf"), pc, device = cairo_pdf,
       width = 17, height = 9)
ggsave(file.path(OUT, "Fig_umap_labelled_by_condition.png"), pc, width = 17,
       height = 9, dpi = 200, bg = "white")

pc2 <- DimPlot(obj, reduction = REDUC, group.by = "condition", pt.size = 0.2) +
  scale_colour_manual(values = COND, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  ggtitle("Condition") + th + AXES
ggsave(file.path(OUT, "Fig_umap_by_condition.pdf"), pc2, device = cairo_pdf,
       width = 9, height = 9)
ggsave(file.path(OUT, "Fig_umap_by_condition.png"), pc2, width = 9, height = 9,
       dpi = 200, bg = "white")

## ---- d. numbered + key ------------------------------------------------------
# key is ID -> name only; nuclei counts and percentages live in
# cluster_label_sizes.csv rather than on the figure
key <- merge(sz, map[, .(cluster, label)], by = "cluster", all.x = TRUE)
key <- key[order(as.numeric(cluster))]
key[, lbl := sprintf("%s  %s", cluster, label)]
pd <- add_labels(
  DimPlot(obj, reduction = REDUC, group.by = CLUSTER_COL, label = FALSE,
          pt.size = 0.2) +
    scale_colour_manual(values = CL_COLS, labels = key$lbl, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    ggtitle("Clusters with cell-type key") + th +
    theme(legend.text = element_text(size = 12, family = FAM)),
  cen_cluster)
ggsave(file.path(OUT, "Fig_umap_numbered_with_key.pdf"), pd, device = cairo_pdf,
       width = 14, height = 9.5)
ggsave(file.path(OUT, "Fig_umap_numbered_with_key.png"), pd, width = 14, height = 9.5,
       dpi = 200, bg = "white")

cat("\nwrote ->", OUT, "\n")
cat("DONE\n")
