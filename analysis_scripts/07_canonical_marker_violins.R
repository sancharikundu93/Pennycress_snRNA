#!/usr/bin/env Rscript
# Step 07 -- canonical marker VIOLINS + dot plot, to accompany the FeaturePlots.
#
# WHY. A FeaturePlot of a canonical marker in snRNA-seq reads as scatter: most nuclei
# are zero for any single gene (see pct_Control/pct_Salt in the step 05 table -- GL2
# is expressed in 0.9% of nuclei, EXP7 in 0.1%), so a UMAP coloured by that gene is
# mostly grey with a few dots. That is a property of the data, not a plotting fault,
# and it is why marker validation is normally shown grouped by cell type rather than
# spatially. These are the readable companions.
#
# THREE FIGURES
#   1. stacked violin  -- markers (rows) x cell types (columns). The standard marker
#                         panel: does each marker peak in the cell type it should?
#   2. dot plot        -- % expressing (size) x scaled mean expression (colour).
#                         Handles sparse genes better than a violin, since a violin
#                         of a 0.9%-expressed gene is a flat line with a spike.
#   3. per-marker violin, split Control vs Salt, one page per marker.
#
# Ortholog choice is the same rule as step 05: a TAIR symbol maps to 1..n MN106 loci
# (family-level for ~77% of genes); the highest-expressed ortholog is plotted and
# n_orthologs is carried into the label so a 9-ortholog "marker" is not mistaken for
# a gene-level result.
#
#   MARKER_SET=literature sbatch run_07_marker_violins.sh   # default
#   MARKER_SET=all        sbatch --export=ALL run_07_marker_violins.sh
#
# Outputs -> ncount500_08132026/results/07_marker_violins/

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
options(bitmapType = "cairo")

FAM  <- "Times New Roman"
COND <- c(Control = "#F0BC35", Salt = "#C2185B")
BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
LAB  <- file.path(ROOT, "results/02_annotation/labels_res0.6.csv")
OUT  <- file.path(ROOT, "results/07_marker_violins")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PANEL <- file.path(BASE, "results/marker_featureplots_res0.6_073126/marker_panel.csv")
ORTH  <- file.path(BASE, "results/marker_featureplots_res0.6_073126/panel_resolved_orthologs.csv")
stopifnot(file.exists(RDS), file.exists(LAB), file.exists(PANEL), file.exists(ORTH))

MARKER_SET  <- Sys.getenv("MARKER_SET", "literature")
CLUSTER_COL <- Sys.getenv("CLUSTER_COL", "RNA_snn_res.0.6")

panel <- fread(PANEL)
orth  <- fread(ORTH); setnames(orth, tolower(names(orth)))
if (MARKER_SET == "literature") panel <- panel[source == "literature"]
cat("markers in set:", nrow(panel), "\n")

obj <- readRDS(RDS)
lab <- fread(LAB, colClasses = "character")
cl  <- as.character(obj@meta.data[[CLUSTER_COL]])
hit <- match(cl, lab$cluster)
if (anyNA(hit)) stop("unlabelled cluster(s): ",
                     paste(sort(unique(cl[is.na(hit)])), collapse = ", "))
# Several clusters share a cell-type name (3 Endodermis, 2 Pericycle, 2 Atrichoblast,
# 2 S phase, 2 Cortex). Keep the cluster ID in the axis label so a marker peaking in
# only ONE of the two Pericycle clusters stays visible instead of being averaged away.
obj$cell_type  <- lab$cell_type[hit]
obj$ct_cluster <- sprintf("%s (%s)", obj$cell_type, cl)
obj$condition  <- factor(as.character(obj$condition), levels = names(COND))
DefaultAssay(obj) <- "RNA"
cat("nuclei:", ncol(obj), "\n")

## ---- resolve symbols to features -------------------------------------------
orth <- orth[feature %in% rownames(obj)][symbol %in% panel$symbol]
expr <- GetAssayData(obj, layer = "data")
me   <- Matrix::rowMeans(expr[unique(orth$feature), , drop = FALSE])
orth[, mean_expr := me[feature]]
setorder(orth, symbol, -mean_expr)
chosen <- orth[, .(feature = feature[1], n_orthologs = .N), by = symbol]
chosen <- merge(panel, chosen, by = "symbol")[!is.na(feature)]
setorder(chosen, cell_type, symbol)
# a symbol whose best ortholog is shared with another symbol is not independent
# evidence (XCP1 and XCP2 both resolve to TaMN106.1G252900) -- say so in the label
chosen[, dup_feature := .N > 1, by = feature]
chosen[, lbl := sprintf("%s (%s)%s", symbol, cell_type,
                        ifelse(n_orthologs > 1, sprintf(" [%d orth]", n_orthologs), ""))]
cat("plotting", nrow(chosen), "markers\n")
if (any(chosen$dup_feature)) {
  cat("NOTE: symbols sharing one MN106 locus (not independent):\n")
  print(chosen[dup_feature == TRUE, .(symbol, feature)])
}
fwrite(chosen[, .(symbol, tair, cell_type, source, feature, n_orthologs, dup_feature)],
       file.path(OUT, "markers_plotted.csv"))

# DEDUPLICATE. Two symbols can resolve to the SAME MN106 locus (XCP1 and XCP2 both
# give TaMN106.1G252900), and Seurat's DotPlot builds a factor over the feature
# vector -- a repeated feature makes scale_x_discrete fail with
# "factor level [n] is duplicated" (job 668261). Collapse them to one entry whose
# label names both symbols, which is also the honest presentation: it is one
# pennycress gene, not two independent markers.
chosen <- chosen[, .(symbol = paste(symbol, collapse = "/"),
                     tair = paste(tair, collapse = "/"),
                     cell_type = cell_type[1], source = source[1],
                     n_orthologs = n_orthologs[1], dup_feature = .N > 1),
                 by = feature]
# Feature names become the row labels below. Seurat silently rewrites underscores to
# dashes in feature names, and cell_type here contains them (Xylem_pole_pericycle),
# so build the label already normalised rather than letting that happen unpredictably.
# Square brackets are dropped for the same reason.
chosen[, lbl := sprintf("%s (%s)%s", symbol, gsub("_", " ", cell_type),
                        ifelse(n_orthologs > 1, sprintf(" x%d", n_orthologs), ""))]
# The SYMBOL can carry an underscore too (WOL_CRE1), and normalising only cell_type
# left those through. CreateAssayObject then rewrites the stored row name to
# "WOL-CRE1 ..." while chosen$lbl still says "WOL_CRE1 ...", so every later lookup by
# label misses: the per-marker loop dies outright ("None of the requested variables
# were found", job 668730) and -- worse, because it is silent -- VlnPlot(stack) and
# DotPlot simply DROP those markers, since Seurat only errors when NO feature matches.
# Apply Seurat's own substitution up front so the stored name and the label agree.
chosen[, lbl := gsub("_", "-", lbl)]
setorder(chosen, cell_type, symbol)
stopifnot(!any(duplicated(chosen$feature)), !any(duplicated(chosen$lbl)))
cat("unique features after collapsing shared loci:", nrow(chosen), "\n")

feats <- chosen$feature
names(feats) <- chosen$lbl

# order cell types by size so the eye starts where the evidence is
ord <- names(sort(table(obj$ct_cluster), decreasing = TRUE))
obj$ct_cluster <- factor(obj$ct_cluster, levels = ord)
Idents(obj) <- "ct_cluster"

th <- theme_classic(base_family = FAM, base_size = 10) +
  theme(text = element_text(colour = "black", family = FAM),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(family = FAM, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(family = FAM, colour = "black"))

## ---- build a marker assay whose ROW NAMES are the symbols -------------------
# VlnPlot/DotPlot label rows by feature name, so plotting the raw MN106 IDs makes an
# unreadable figure (job 668261/668262 produced exactly that: strips reading
# "TaMN106.3G288300.v4.1"). Patching strip text after the fact is fragile because
# stack=TRUE builds facets internally. Instead give Seurat a small assay that is
# already named the way the figure should read.
mat <- expr[chosen$feature, , drop = FALSE]
rownames(mat) <- chosen$lbl
obj[["MARKERS"]] <- CreateAssayObject(data = mat)
DefaultAssay(obj) <- "MARKERS"
plot_feats <- chosen$lbl
# Fail loudly rather than let a figure ship with markers silently missing.
missing <- setdiff(plot_feats, rownames(obj[["MARKERS"]]))
if (length(missing)) stop("label(s) not stored verbatim in the MARKERS assay: ",
                          paste(missing, collapse = ", "))

## ---- 1. stacked violin ------------------------------------------------------
# VlnPlot(stack=TRUE) puts markers on rows sharing one x axis -- the compact form
# that fits on a slide, unlike 45 separate violins.
p1 <- VlnPlot(obj, features = plot_feats, stack = TRUE, flip = TRUE,
              fill.by = "feature", pt.size = 0) +
  NoLegend() + th +
  labs(title = "Canonical markers by cell type", x = NULL, y = NULL)
ggsave(file.path(OUT, paste0("Fig_violin_stacked_", MARKER_SET, ".pdf")), p1,
       device = cairo_pdf, width = 11, height = 1 + 0.55 * length(feats),
       limitsize = FALSE)
ggsave(file.path(OUT, paste0("Fig_violin_stacked_", MARKER_SET, ".png")), p1,
       width = 11, height = 1 + 0.55 * length(feats), dpi = 200, bg = "white",
       limitsize = FALSE)

## ---- 2. dot plot ------------------------------------------------------------
p2 <- DotPlot(obj, features = plot_feats, dot.scale = 5) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                         midpoint = 0, name = "scaled\nexpression") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Canonical markers: % expressing and scaled mean expression",
       x = NULL, y = NULL) + coord_flip()
ggsave(file.path(OUT, paste0("Fig_dotplot_", MARKER_SET, ".pdf")), p2,
       device = cairo_pdf, width = 10, height = 1 + 0.32 * length(feats),
       limitsize = FALSE)
ggsave(file.path(OUT, paste0("Fig_dotplot_", MARKER_SET, ".png")), p2,
       width = 10, height = 1 + 0.32 * length(feats), dpi = 200, bg = "white",
       limitsize = FALSE)

## ---- 3. per-marker violin, Control vs Salt ----------------------------------
mk_one <- function(i) {
  VlnPlot(obj, features = chosen$lbl[i], split.by = "condition",
          pt.size = 0, cols = unname(COND)) +
    scale_fill_manual(values = COND, name = NULL) +
    th + labs(title = chosen$lbl[i], x = NULL, y = "expression")
}
gs <- lapply(seq_len(nrow(chosen)), mk_one)
g3 <- wrap_plots(gs, ncol = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "top")
ggsave(file.path(OUT, paste0("Fig_violin_by_condition_", MARKER_SET, ".pdf")), g3,
       device = cairo_pdf, width = 11, height = 2.6 * nrow(chosen), limitsize = FALSE)

cat("\nwrote ->", OUT, "\n")
cat("DONE\n")
