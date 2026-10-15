#!/usr/bin/env Rscript
# =============================================================================
# Step 75 -- labelled UMAP on the step-74 from-scratch annotation
#
#   Rscript 75_umap_relabelled.R
#
# Renders the res 0.6 UMAP carrying the labels re-derived in step 74 from THIS object,
# with no inheritance from the 1000-gate annotation. Two labels change:
#   cluster 15  "Seed/embryo remnant"  -> "LEA/desiccation-high"
#   cluster 22  "Photosynthetic cell"  -> "Photosynthetic epidermis"
# Everything else is unchanged, so this figure is directly comparable to
# 06_labelled_umap/Fig_umap_labelled.png.
#
# NO OBJECT LOAD. Coordinates come from 61_auxin_priming_umap's saved UMAP table, which
# carries the cell barcode, so the join to 51_doublets (barcode -> cluster) is BY KEY and
# not positional -- unlike steps 73/74, no order assumption is needed or made here. The
# script asserts a complete join rather than silently dropping unmatched nuclei.
#
# Outputs -> ncount500_08132026/results/75_umap_relabelled/
# =============================================================================
.libPaths(c("/home/sancharikundu/R/tinkercliffs-rome/4.4.2",
            "/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
options(bitmapType = "cairo")

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
LIN  <- file.path(BASE, "ncount500_08132026")
XY   <- file.path(LIN, "results/61_auxin_priming_umap/auxin_LR_score_umap_coords.csv.gz")
DBL  <- file.path(LIN, "results/51_doublets/per_nucleus_doublet_calls.csv.gz")
NEW  <- file.path(LIN, "results/74_annotation_from_scratch/labels_res0.6_from_scratch.csv")
OUT  <- file.path(LIN, "results/75_umap_relabelled")
source(file.path(LIN, "scripts/palette_celltypes.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FAM <- "Times New Roman"
cat("=== step 75: relabelled UMAP ===\n")
for (f in c(XY, DBL, NEW)) stopifnot(file.exists(f))

xy  <- fread(XY)[, .(cell, x, y, condition)]
dbl <- fread(DBL)[, .(cell, cluster)]
new <- fread(NEW)[, .(cluster = as.integer(cluster), cell_type, confidence)]

d <- merge(xy, dbl, by = "cell")
stopifnot(nrow(d) == nrow(xy))          # a keyed join must not lose nuclei
d <- merge(d, new, by = "cluster")
stopifnot(nrow(d) == nrow(xy), !anyNA(d$cell_type))
cat("joined", nrow(d), "nuclei by barcode; no losses\n")

## The two renamed clusters inherit their predecessors' palette entries, so the figure is
## colour-comparable to every earlier UMAP -- only the legend text changes.
PAL2 <- PAL
PAL2["LEA/desiccation-high"]     <- unname(PAL["Seed/embryo remnant"])
PAL2["Photosynthetic epidermis"] <- unname(PAL["Photosynthetic cell"])
LV <- c(CT_LEVELS[CT_LEVELS %in% unique(d$cell_type)],
        intersect(c("LEA/desiccation-high", "Photosynthetic epidermis"), unique(d$cell_type)))
stopifnot(setequal(LV, unique(d$cell_type)))
d[, cell_type := factor(cell_type, levels = LV)]

## orient_xy returns an UNNAMED list (palette_celltypes.R), so index by position --
## $x/$y silently yields NULL and data.table then tries to DELETE the column.
## Step 61 saved raw coordinates (it never calls orient_xy), so the shared ORIENT is
## applied here, which is what keeps this UMAP in the same orientation as every other
## figure in the folder.
xyo <- orient_xy(d$x, d$y)
d[, X := xyo[[1]]][, Y := xyo[[2]]]
lab <- d[, .(X = median(X), Y = median(Y), n = .N), by = cell_type][order(cell_type)]
## Several cell types span MORE THAN ONE cluster (Endodermis 3, 11, 17; Atrichoblast 4,
## 13; S phase 7, 10), so joining the per-CLUSTER confidence straight onto a per-CELL-TYPE
## centroid table multiplies the rows and draws the same label two or three times on top
## of itself. Collapse to one row per cell type first, keeping the clusters it covers.
conf <- new[, .(clusters = paste(sort(cluster), collapse = "/"),
                confidence = paste(unique(confidence), collapse = "/")), by = cell_type]
lab <- merge(lab, conf, by = "cell_type")
fwrite(lab, file.path(OUT, "label_centroids.csv"))
cat("\ncell types (", nrow(lab), "):\n"); print(lab[, .(cell_type, clusters, n, confidence)])

base <- function(p, title, sub)
  p + scale_colour_manual(values = PAL2[LV], drop = FALSE, name = NULL) +
    coord_equal() +
    theme_void(base_family = FAM) +
    theme(text = element_text(family = FAM, colour = "black"),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          legend.text = element_text(size = 11)) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    labs(title = title, subtitle = sub)

SUB <- paste0(format(nrow(d), big.mark = ","),
              " nuclei, resolution 0.6, nCount 500-30000. Labels re-derived from this ",
              "object (step 74) — no transfer from the 1000-gate annotation.")

p1 <- base(ggplot(d, aes(X, Y, colour = cell_type)) +
             geom_point(size = 0.25, alpha = 0.75, show.legend = TRUE) +
             ggrepel::geom_text_repel(data = lab, aes(X, Y, label = cell_type),
                                      inherit.aes = FALSE, family = FAM, size = 4.2,
                                      colour = "black", bg.colour = "white", bg.r = 0.14,
                                      box.padding = 0.4, max.overlaps = Inf, seed = 1),
           "Pennycress root snRNA-seq atlas", SUB) + theme(legend.position = "none")
ggsave(file.path(OUT, "Fig_umap_relabelled.pdf"), p1, width = 10, height = 9.5,
       device = cairo_pdf)
ggsave(file.path(OUT, "Fig_umap_relabelled.png"), p1, width = 10, height = 9.5,
       dpi = 300, bg = "white")

p2 <- base(ggplot(d, aes(X, Y, colour = cell_type)) + geom_point(size = 0.25, alpha = 0.75),
           "Pennycress root snRNA-seq atlas", SUB)
ggsave(file.path(OUT, "Fig_umap_relabelled_legend.pdf"), p2, width = 12.5, height = 9.5,
       device = cairo_pdf)
ggsave(file.path(OUT, "Fig_umap_relabelled_legend.png"), p2, width = 12.5, height = 9.5,
       dpi = 300, bg = "white")

p3 <- base(ggplot(d, aes(X, Y, colour = cell_type)) + geom_point(size = 0.22, alpha = 0.75) +
             facet_wrap(~ condition),
           "Control vs Salt", SUB) +
  theme(strip.text = element_text(face = "bold", size = 14, family = FAM))
ggsave(file.path(OUT, "Fig_umap_relabelled_by_condition.pdf"), p3, width = 16, height = 8.5,
       device = cairo_pdf)
ggsave(file.path(OUT, "Fig_umap_relabelled_by_condition.png"), p3, width = 16, height = 8.5,
       dpi = 300, bg = "white")

cat("\nwrote 3 figures ->", OUT, "\n")
cat("done:", format(Sys.time()), "\n")
