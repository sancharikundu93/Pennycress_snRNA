#!/usr/bin/env Rscript
# =============================================================================
# Step 78 -- canonical marker FEATURE PLOT grid, thresholded
#
#   sbatch run_78_featureplot_grid.sh
#
# REPLACES 05_marker_featureplots/Fig_canonical_markers_literature.png, which is unusable:
# 1275 x 5864 px (a single narrow column of 17 stacked rows), drawn in the UN-rotated UMAP
# orientation that no other figure in this folder uses, with locus ids printed sideways down
# the margin and two unexplained columns per gene (the two MN106 loci best-hitting one
# Arabidopsis marker). At print size each panel is a few hundred pixels wide.
#
# WHAT CHANGES, BEYOND LAYOUT.
# A raw feature plot colours every nucleus with a non-zero count, so a gene detected at 1-2
# molecules in half the atlas looks "expressed everywhere" and a genuinely specific gene
# looks similar. This one shows only where a marker is STRONG: expression is scaled 0-1 per
# gene against that gene's own 99th percentile (not its max -- one outlier nucleus would
# otherwise compress every other value toward zero), and nuclei below MIN_SCALED are drawn
# in the same pale grey as the non-expressing ones.
#
# THE FILTER DOES NOT LAUNDER A BAD MARKER, and that is the property that matters. Scaling
# is per GENE, so a gene strong in one cluster keeps one coloured island while a gene strong
# in fifteen keeps fifteen. BGLU23 is deliberately included: it is the atrichoblast marker
# and it stays visibly smeared across the whole atlas, exactly as it should.
#
# COORDINATES come from step 61's saved UMAP table, joined to the object BY BARCODE, and the
# shared ORIENT is applied so this figure matches every other UMAP in the folder. The join is
# asserted complete rather than left to silently drop nuclei.
#
# Outputs -> ncount500_08132026/results/78_featureplot_grid/
# =============================================================================
.libPaths(c("/home/sancharikundu/R/tinkercliffs-rome/4.4.2",
            "/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); library(Matrix); library(patchwork)
})
options(bitmapType = "cairo")

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
LIN  <- file.path(BASE, "ncount500_08132026")
RES  <- "0.6"
RDS  <- file.path(LIN, "results/penny_mn106_ncount500_integrated.rds")
XY   <- file.path(LIN, "results/61_auxin_priming_umap/auxin_LR_score_umap_coords.csv.gz")
DBL  <- file.path(LIN, "results/51_doublets/per_nucleus_doublet_calls.csv.gz")
## The symbol -> MN106 locus map. 05_marker_featureplots/marker_table_for_reference.csv
## holds only 17 symbols (the subset step 05 happened to draw); the FULL 44-marker panel
## resolves through the same two files step 02 uses, so the gene list here is not limited
## to whatever the old figure chose. Where a symbol pulls several loci the most-expressed
## one is taken, matching step 02's rule, and n_orthologs is carried so a family-level
## "marker" is visible as such.
## NOTE the res0.8 path: results/marker_featureplots_res0.6_073126/ does not exist on disk,
## only the res0.8 sibling does. That is fine and not a resolution mismatch -- this file is
## a symbol -> MN106 locus map (32 symbols, 140 rows) with no cluster information in it, so
## it is resolution-independent. Step 02 points at the res0.6 path and silently skips the
## whole canonical block when it is absent, which is worth knowing.
ORTH  <- file.path(BASE,
           "results/marker_featureplots_res0.8_073126/panel_resolved_orthologs.csv")
LAB  <- Sys.getenv("LABELS", file.path(LIN, "results/02_annotation/labels_res0.6.csv"))
OUT  <- file.path(LIN, "results/78_featureplot_grid")
source(file.path(LIN, "scripts/palette_celltypes.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FAM <- "Times New Roman"; INK <- "black"
OFF <- "grey88"                    # below threshold / not expressing
ON  <- c("#DCD3EC", "#4A1486")     # low -> high, above threshold
MIN_SCALED <- as.numeric(Sys.getenv("MIN_SCALED", "0.5"))
DROP <- as.integer(strsplit(Sys.getenv("DROP_CLUSTERS", "5,15"), ",")[[1]])

## Curated so every major tissue gets a panel, ordered anatomically. BGLU23 and NPF7.3 are
## in on purpose as the negative controls -- a reader should SEE what a non-specific marker
## looks like next to CASP1, or the panel is just an assertion.
GENES <- strsplit(Sys.getenv("GENES", paste(
  "WOX5,ACR4,PRP3,AGP3,EXP7,COBL9,BGLU23,GL2,ROSY1,CASP1,ESB1,SCR,LTPG5,",
  "IQD13,NPF7.3,SUC2,APL,XCP1,XCP2", sep = "")), ",")[[1]]
GENES <- trimws(GENES[nzchar(trimws(GENES))])

cat("=== step 78: thresholded feature plot grid ===\n")
cat("genes:", paste(GENES, collapse = ", "), "\n")
cat("MIN_SCALED =", MIN_SCALED, "\n")
for (f in c(RDS, XY, DBL, ORTH, LAB)) stopifnot(file.exists(f))

orth <- fread(ORTH); setnames(orth, tolower(names(orth)))
miss <- setdiff(GENES, orth$symbol)
## Warn and continue rather than stop: an unavailable symbol should cost one panel, not the
## whole figure. The names are printed so the omission is never silent.
if (length(miss)) {
  cat("NOT resolvable to an MN106 locus, dropped:", paste(miss, collapse = ", "), "\n")
  cat("available symbols:", paste(sort(unique(orth$symbol)), collapse = ", "), "\n")
  GENES <- setdiff(GENES, miss)
}
stopifnot(length(GENES) > 0)

xy  <- fread(XY)[, .(cell, x, y)]
dbl <- fread(DBL)[, .(cell, cluster)]
lab <- fread(LAB)[, .(cluster = as.integer(cluster), cell_type)]
md  <- merge(merge(xy, dbl, by = "cell"), lab, by = "cluster")
stopifnot(nrow(md) == nrow(xy))
md <- md[!cluster %in% DROP]
cat("nuclei after dropping clusters", paste(DROP, collapse = "/"), ":", nrow(md), "\n")

xyo <- orient_xy(md$x, md$y); md[, X := xyo[[1]]][, Y := xyo[[2]]]

cat("loading object ...\n")
so <- readRDS(RDS)
ex  <- GetAssayData(so, assay = "RNA", layer = "data")
o   <- orth[symbol %in% GENES][feature %in% rownames(ex)]
if (!length(unique(o$symbol))) stop("no requested symbol has a locus in the matrix")
o[, mean_expr := Matrix::rowMeans(ex[feature, , drop = FALSE])[feature]]
setorder(o, symbol, -mean_expr)
chosen <- o[, .(feature = feature[1], n_orthologs = .N), by = symbol]
gone <- setdiff(GENES, chosen$symbol)
if (length(gone)) { cat("dropped, no locus in the matrix:", paste(gone, collapse = ", "), "\n")
                    GENES <- setdiff(GENES, gone) }
g2f <- setNames(chosen$feature, chosen$symbol)[GENES]
cat("\nresolved loci:\n"); print(chosen[match(GENES, symbol)])
E <- ex[unname(g2f), , drop = FALSE]
rm(ex)
E <- E[, md$cell, drop = FALSE]          # reorder to the plotting table, by barcode
stopifnot(identical(colnames(E), md$cell))
rm(so); invisible(gc())

## Reference UMAP, same palette and orientation as every other figure.
ct_lv <- CT_LEVELS[CT_LEVELS %in% unique(md$cell_type)]
md[, cell_type := factor(cell_type, levels = ct_lv)]
cen <- md[, .(X = median(X), Y = median(Y)), by = cell_type]
p_ref <- ggplot(md, aes(X, Y, colour = cell_type)) +
  geom_point(size = 0.18, alpha = 0.7, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = cen, aes(X, Y, label = cell_type), inherit.aes = FALSE,
                           family = FAM, size = 3.1, colour = INK, bg.colour = "white",
                           bg.r = 0.13, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = PAL[ct_lv]) + coord_equal() + theme_void(base_family = FAM) +
  labs(title = "Cell type") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = FAM))

summ <- list()
panels <- lapply(seq_along(GENES), function(i) {
  v <- as.numeric(E[g2f[GENES[i]], ])
  q <- stats::quantile(v, 0.99, names = FALSE)
  sc <- if (q > 0) pmin(v / q, 1) else rep(0, length(v))
  keep <- sc >= MIN_SCALED
  summ[[i]] <<- data.table(gene = GENES[i], feature = g2f[GENES[i]],
                           q99 = round(q, 3), n_shown = sum(keep),
                           pct_shown = round(100 * mean(keep), 2))
  d <- copy(md)[, `:=`(val = fifelse(keep, sc, NA_real_))]
  setorder(d, val, na.last = FALSE)      # coloured nuclei drawn on top of the grey ones
  ggplot(d, aes(X, Y)) +
    geom_point(data = d[is.na(val)], colour = OFF, size = 0.16, alpha = 0.55) +
    geom_point(data = d[!is.na(val)], aes(colour = val), size = 0.34) +
    scale_colour_gradient(low = ON[1], high = ON[2], limits = c(MIN_SCALED, 1),
                          name = "Scaled\nexpression") +
    coord_equal() + theme_void(base_family = FAM) +
    labs(title = GENES[i]) +
    theme(plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5,
                                    family = FAM),
          legend.position = "none")
})

sm <- rbindlist(summ)
fwrite(sm, file.path(OUT, "featureplot_summary.csv"))
cat("\n=== nuclei shown per gene (above threshold) ===\n"); print(sm)

NC <- as.integer(Sys.getenv("NCOL", "5"))
g <- wrap_plots(c(list(p_ref), panels), ncol = NC) +
  plot_annotation(
    title = "Canonical marker expression across the root atlas",
    subtitle = paste0("Expression scaled 0-1 per gene against its own 99th percentile; ",
                      "nuclei below ", MIN_SCALED, " are grey.\nBGLU23 and NPF7.3 are ",
                      "included as non-specific controls, COBL9 as a failed one. Clusters ",
                      paste(DROP, collapse = " and "), " are not shown."),
    theme = theme(plot.title = element_text(family = FAM, face = "bold", size = 17),
                  plot.subtitle = element_text(family = FAM, size = 11)))

NR <- ceiling((length(GENES) + 1) / NC)
W <- 3.0 * NC; H <- 3.0 * NR + 1.0
f <- file.path(OUT, "Fig_marker_featureplots")
ggsave(paste0(f, ".pdf"), g, width = W, height = H, device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(f, ".png"), g, width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("\nwrote %s (%.1f x %.1f in, %d panels)\n", f, W, H, length(GENES) + 1))
cat("done:", format(Sys.time()), "\n")
