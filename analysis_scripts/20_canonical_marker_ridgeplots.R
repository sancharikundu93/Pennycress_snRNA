#!/usr/bin/env Rscript
# Step 20 -- canonical marker RIDGE PLOTS, for the supplementary marker-validation
# figure. Completes the set: step 05 FeaturePlots (spatial), step 07 violins + dot
# plot (grouped), step 20 ridges (per-cell-type expression distributions).
#
# WHY A RIDGE ADDS ANYTHING OVER THE VIOLIN. A violin is symmetric about its axis, so
# a marker expressed in a minority of nuclei reads as a thin spike either way and the
# eye cannot judge how much of the cell type carries it. A ridge is a one-sided
# density on a shared x axis, so the cell types stack and can be compared directly --
# which is the claim a marker panel has to support: this gene's mass sits in THIS
# cell type and not the ones above and below it.
#
# RIDGES ARE BUILT FROM EXPRESSING NUCLEI ONLY (RIDGE_MODE=expressing, default).
# The first run of this step plotted all nuclei and produced 15 panels of flat lines:
# most canonical markers are detected in a small percentage of nuclei (GL2 0.9%,
# EXP7 0.1% -- step 05's table), so the density is a point mass at zero, the
# bandwidth collapses to ~1e-06, and the real expression is invisible. Conditioning
# on detection is the standard fix and asks the question the panel is actually for:
# WHERE a marker is detected, how strongly is it expressed?
#
# The cost is that a ridge no longer says how MANY nuclei carry the gene, so the
# percentage expressing is printed at the right of every ridge and a cell type with
# fewer than MIN_N expressing nuclei is left blank rather than drawn from noise.
# Read the two together: a tall ridge at 0.4% is a handful of nuclei.
# RIDGE_MODE=all reproduces the all-nuclei version if it is ever wanted.
#
# The marker set, the ortholog choice and the label construction are IDENTICAL to
# step 07 -- same panel file, same "highest-expressed ortholog" rule, same collapsing
# of symbols that share one MN106 locus -- so the three figures show the same genes
# under the same names and can be cited as one panel.
#
#   MARKER_SET=literature sbatch run_20_marker_ridgeplots.sh   # default
#   MARKER_SET=all        sbatch --export=ALL run_20_marker_ridgeplots.sh
#
# Outputs -> ncount500_08132026/results/20_marker_ridgeplots/
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
  library(ggridges)
})
options(bitmapType = "cairo")

FAM  <- "Times New Roman"
BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
LAB  <- file.path(ROOT, "results/02_annotation/labels_res0.6.csv")
OUT  <- file.path(ROOT, "results/20_marker_ridgeplots")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/palette_celltypes.R"))

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
# Keep the cluster ID on the axis, as step 07 does: five cell-type names are carried
# by more than one cluster (3 Endodermis, 2 Atrichoblast, 2 S phase, 2 Cortex), and a
# marker that peaks in only one of them is exactly what this figure has to show.
obj$cell_type  <- lab$cell_type[hit]
obj$ct_cluster <- sprintf("%s (%s)", obj$cell_type, cl)
DefaultAssay(obj) <- "RNA"
cat("nuclei:", ncol(obj), "\n")

## ---- resolve symbols to features (same rule as step 07) ---------------------
orth <- orth[feature %in% rownames(obj)][symbol %in% panel$symbol]
expr <- GetAssayData(obj, layer = "data")
me   <- Matrix::rowMeans(expr[unique(orth$feature), , drop = FALSE])
orth[, mean_expr := me[feature]]
setorder(orth, symbol, -mean_expr)
chosen <- orth[, .(feature = feature[1], n_orthologs = .N), by = symbol]
chosen <- merge(panel, chosen, by = "symbol")[!is.na(feature)]
setorder(chosen, cell_type, symbol)

# Two symbols can resolve to the SAME MN106 locus (XCP1/XCP2 -> TaMN106.1G252900).
# Collapse them: it is one pennycress gene, not two independent markers, and a
# duplicated feature also breaks the discrete scale downstream.
chosen <- chosen[, .(symbol = paste(symbol, collapse = "/"),
                     tair = paste(tair, collapse = "/"),
                     cell_type = cell_type[1], source = source[1],
                     n_orthologs = n_orthologs[1], dup_feature = .N > 1),
                 by = feature]
# Seurat rewrites underscores in feature names to dashes when an assay is created, so
# normalise the label up front -- otherwise every later lookup by label silently
# misses and markers vanish from the figure without an error (step 07's job 668730).
chosen[, lbl := sprintf("%s (%s)%s", symbol, gsub("_", " ", cell_type),
                        ifelse(n_orthologs > 1, sprintf(" x%d", n_orthologs), ""))]
chosen[, lbl := gsub("_", "-", lbl)]
setorder(chosen, cell_type, symbol)
stopifnot(!any(duplicated(chosen$feature)), !any(duplicated(chosen$lbl)))
cat("unique features after collapsing shared loci:", nrow(chosen), "\n")
fwrite(chosen[, .(symbol, tair, cell_type, source, feature, n_orthologs, dup_feature)],
       file.path(OUT, paste0("markers_plotted_", MARKER_SET, ".csv")))

## ---- expression matrix, rows named as the panels should read ----------------
# Panels are titled from these row names; the raw IDs would title them
# "TaMN106.3G288300.v4.1".
mat <- expr[chosen$feature, , drop = FALSE]
rownames(mat) <- chosen$lbl

## ---- identity order: anatomical, so every panel reads the same way ----------
# Ridges stack bottom-to-top, so ordering by CT_LEVELS makes each panel read outer
# tissue -> stele exactly like the proportion bar, instead of by cluster number.
ct_order <- unique(obj$ct_cluster)[
  order(match(sub(" \\(.*", "", unique(obj$ct_cluster)), CT_LEVELS),
        as.integer(sub(".*\\((\\d+)\\)", "\\1", unique(obj$ct_cluster))))]
obj$ct_cluster <- factor(obj$ct_cluster, levels = rev(ct_order))
Idents(obj) <- "ct_cluster"
n_id <- nlevels(obj$ct_cluster)
cat("identities:", n_id, "\n")

# One colour per CELL TYPE from the project palette, so a ridge matches its cluster
# in the UMAP. Clusters sharing a cell type deliberately share a colour -- the axis
# label carries the cluster ID that separates them.
id_cols <- ct_palette(sub(" \\(.*", "", levels(obj$ct_cluster)))[
  sub(" \\(.*", "", levels(obj$ct_cluster))]
names(id_cols) <- levels(obj$ct_cluster)

th <- theme_classic(base_family = FAM, base_size = 10) +
  theme(text = element_text(colour = "black", family = FAM),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(family = FAM, face = "bold", size = 9),
        plot.margin = margin(4, 34, 4, 4),
        legend.position = "none")

RIDGE_MODE <- Sys.getenv("RIDGE_MODE", "expressing")
if (!RIDGE_MODE %in% c("expressing", "all"))
  stop("RIDGE_MODE must be 'expressing' or 'all', got: ", RIDGE_MODE)
# Below this many detected nuclei a density is noise dressed as a distribution, so
# the row is left blank -- the percentage beside it still reports what was there.
MIN_N <- as.integer(Sys.getenv("MIN_N", "10"))

ids <- obj$ct_cluster   # factor, levels already in plotting order
stopifnot(identical(colnames(mat), colnames(obj)))

# Detection rate per cell type per marker, written out alongside the figure: it is
# the number the ridges cannot show, and the one a reviewer will ask for.
det <- rbindlist(lapply(chosen$lbl, function(l) {
  v <- as.numeric(mat[l, ])
  data.table(marker = l, id = ids, x = v)[
    , .(n = .N, n_expr = sum(x > 0), pct_expr = 100 * mean(x > 0),
        mean_expr_in_expressing = mean(x[x > 0])), by = id][, marker := l][]
}))
fwrite(det[, .(marker, cell_type_cluster = id, n, n_expr,
               pct_expr = round(pct_expr, 2),
               mean_expr_in_expressing = round(mean_expr_in_expressing, 3))],
       file.path(OUT, paste0("detection_by_celltype_", MARKER_SET, ".csv")))

mk_ridge <- function(lbl) {
  d <- data.table(id = ids, x = as.numeric(mat[lbl, ]))
  st <- d[, .(pct = 100 * mean(x > 0), n_expr = sum(x > 0)), by = id]
  dd <- if (RIDGE_MODE == "expressing") d[x > 0] else d
  dd <- dd[id %in% st[n_expr >= MIN_N, id]]
  xmax <- if (nrow(dd)) max(dd$x) else 1
  p <- ggplot(dd, aes(x = x, y = id, fill = id)) +
    geom_density_ridges(scale = 1.7, rel_min_height = 0.008, alpha = 0.9,
                        colour = "grey30", linewidth = 0.18) +
    # every cell type keeps its row even with no ridge, so panels stay comparable
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(values = id_cols, drop = FALSE) +
    # percentage detected, printed just outside the panel
    geom_text(data = st, aes(x = xmax * 1.06, y = id,
                             label = sprintf("%.1f%%", pct)),
              inherit.aes = FALSE, hjust = 0, size = 2.1, family = FAM,
              colour = "grey25") +
    coord_cartesian(xlim = c(0, xmax * 1.02), clip = "off") +
    th + labs(title = lbl, x = "expression (detected nuclei)", y = NULL)
  if (RIDGE_MODE == "all") p <- p + labs(x = "expression (all nuclei)")
  p
}

## ---- 1. one grid over every marker ------------------------------------------
NC  <- as.integer(Sys.getenv("NCOL", "3"))
SFX <- if (RIDGE_MODE == "all") paste0(MARKER_SET, "_allnuclei") else MARKER_SET
gs <- lapply(chosen$lbl, mk_ridge)
g  <- wrap_plots(gs, ncol = NC)
nr <- ceiling(length(gs) / NC)
W  <- 4.6 * NC
H  <- (0.9 + 0.20 * n_id) * nr
ggsave(file.path(OUT, paste0("Fig_ridge_", SFX, ".pdf")), g,
       device = cairo_pdf, width = W, height = H, limitsize = FALSE)
ggsave(file.path(OUT, paste0("Fig_ridge_", SFX, ".png")), g,
       width = W, height = H, dpi = 200, bg = "white", limitsize = FALSE)
cat("wrote grid:", length(gs), "markers,", nr, "rows x", NC, "cols\n")

## ---- 2. one file per marker cell type ---------------------------------------
# Same structure as step 05's per-cell-type FeaturePlot files, so a supplementary
# figure can cite the endodermis markers without carrying all 45 panels.
for (ct in unique(chosen$cell_type)) {
  sub <- chosen[cell_type == ct]
  gg  <- wrap_plots(lapply(sub$lbl, mk_ridge), ncol = min(2, nrow(sub)))
  nrs <- ceiling(nrow(sub) / min(2, nrow(sub)))
  ggsave(file.path(OUT, sprintf("Fig_ridge_%s_%s.pdf", gsub("[^A-Za-z0-9]+", "_", ct),
                                SFX)),
         gg, device = cairo_pdf, width = 4.6 * min(2, nrow(sub)),
         height = (0.9 + 0.20 * n_id) * nrs, limitsize = FALSE)
  cat("  ", ct, ":", nrow(sub), "markers\n")
}

cat("\nwrote ->", OUT, "\n")
cat("DONE\n")
