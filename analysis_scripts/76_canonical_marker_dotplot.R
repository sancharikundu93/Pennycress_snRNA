#!/usr/bin/env Rscript
# =============================================================================
# Step 76 -- canonical marker DOT PLOT: markers x clusters, the annotation-justification
#            figure
#
#   Rscript 76_canonical_marker_dotplot.R
#
# WHY A DOT PLOT AND NOT MORE FEATURE PLOTS.
# 05_marker_featureplots/Fig_canonical_markers_literature.png is 1275 x 5864 px -- a single
# narrow column of 17 stacked rows, in the UN-rotated UMAP orientation that no other figure
# in this folder uses, with the locus IDs printed sideways down the margin. Beyond being
# unusable at print size, a feature plot shows WHERE a gene is detected but not HOW
# SPECIFICALLY, and specificity is the entire claim an annotation figure has to make. A
# marker detected in 92% of pericycle nuclei looks convincing on a UMAP until you notice it
# is also >46% in 22 of the other 23 clusters (NPF7.3, and it is not alone).
#
# This panel makes that visible in one glance: dot SIZE = % of nuclei expressing,
# dot COLOUR = mean expression scaled 0-1 WITHIN each marker row. Row-scaling is what makes
# rows comparable -- without it the plot ranks markers by absolute abundance, which is the
# same mistake that broke canonical_best (step 74 header).
#
# NO OBJECT LOAD. 02_annotation/canonical_marker_pct_RNA_snn_res_0_6.csv already stores
# pct_expressing and mean_expr per marker per cluster, which is exactly dot-plot input.
#
# COLUMNS ARE CLUSTERS, NOT CELL TYPES, and that is deliberate. The figure exists to justify
# "cluster 17 IS endodermis", so the cluster id has to be on the axis; collapsing to cell
# type would hide that Endodermis is three clusters (3, 11, 17) which do NOT agree -- CASP1
# is 62.2% in cluster 17 and 4.3% at its next best. Labels read "17 Endodermis".
#
# n_orthologs IS PRINTED ON EVERY ROW, because it predicts failure. Pathway/panel membership
# is Phytozome best-BLAST-hit, not orthology, so a marker pulling many MN106 loci is
# averaging over a gene family. COBL9 (n=5) is the worst case in this panel: a textbook root
# hair marker that lands at 1.0% in Trichoblast -- dead LAST of 23 clusters -- while reaching
# 15.2% in G2/M phase. AGP3 (n=1) hits 55.4% in Trichoblast against 4.7% next best. Read the
# n column before trusting a row.
#
# Outputs -> ncount500_08132026/results/76_marker_dotplot/
# =============================================================================
.libPaths(c("/home/sancharikundu/R/tinkercliffs-rome/4.4.2",
            "/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
options(bitmapType = "cairo")

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
LIN  <- file.path(BASE, "ncount500_08132026")
PCT  <- file.path(LIN, "results/02_annotation/canonical_marker_pct_RNA_snn_res_0_6.csv")
## Original labels (2026-08-20, requested): "Photosynthetic cell" / "Seed/embryo remnant"
## stay for now. Step 74's from-scratch table still disagrees with both and is untouched.
LAB  <- Sys.getenv("LABELS", file.path(LIN, "results/02_annotation/labels_res0.6.csv"))
## Clusters 5 (pericycle_like) and 15 (Seed/embryo remnant) are held out -- both identities
## are still open, so they do not belong in the panel that justifies the annotation.
DROP <- as.integer(strsplit(Sys.getenv("DROP_CLUSTERS", "5,15"), ",")[[1]])
## MIN_SCALED (2026-08-20, requested): hide any dot whose row-scaled mean expression is
## below this, so the panel shows where a marker is STRONG rather than everywhere it is
## detected at all. 0.5 = the midpoint of the colour scale.
##
## This is a DISPLAY filter and it does not launder a bad marker. Scaling is per ROW, so a
## gene high in one cluster keeps one dot while a gene high in fifteen keeps fifteen --
## BGLU23 and NPF7.3 stay visibly non-specific, which is the property the figure must not
## hide. What it removes is the low-level haze that makes every row look occupied.
## Set MIN_SCALED=0 to see every dot again; the unfiltered values are in the CSV either way.
MIN_SCALED <- as.numeric(Sys.getenv("MIN_SCALED", "0.5"))
OUT  <- file.path(LIN, "results/76_marker_dotplot")
source(file.path(LIN, "scripts/palette_celltypes.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FAM <- "Times New Roman"; INK <- "black"
LOW <- "#F2F0F7"; HIGH <- "#4A1486"

cat("=== step 76: canonical marker dot plot ===\n")
stopifnot(file.exists(PCT), file.exists(LAB))
p   <- fread(PCT)
lab <- fread(LAB)[, .(cluster = as.integer(cluster), cell_type)][!cluster %in% DROP]
cat("dropped clusters:", paste(DROP, collapse = ", "), "\n")
d   <- merge(p, lab, by = "cluster")
cat(nrow(d), "marker x cluster rows;", uniqueN(d$symbol), "markers,",
    uniqueN(d$cluster), "clusters\n")

## Columns: clusters ordered anatomically by their cell type, then by cluster id.
ct_ord <- CT_LEVELS[CT_LEVELS %in% unique(lab$cell_type)]
stopifnot(setequal(ct_ord, unique(lab$cell_type)))
lab[, ct_rank := match(cell_type, ct_ord)]
setorder(lab, ct_rank, cluster)
lab[, col_lab := sprintf("%d  %s", cluster, cell_type)]
d <- merge(d, lab[, .(cluster, col_lab, ct_rank)], by = "cluster")
d[, col_lab := factor(col_lab, levels = lab$col_lab)]

## Rows: scale mean_expr to 0-1 within each marker, and record specificity = best / 2nd best
## pct. That ratio is the number the figure is really about, so it goes on the row label.
d[, expr_scaled := {
  r <- range(mean_expr, na.rm = TRUE)
  if (diff(r) == 0) rep(0, .N) else (mean_expr - r[1]) / diff(r)
}, by = symbol]

spec <- d[, {
  s <- sort(pct_expressing, decreasing = TRUE)
  .(best = s[1], second = s[2], ratio = round(s[1] / max(s[2], 0.1), 1),
    best_ct = col_lab[which.max(pct_expressing)],
    panel = marker_cell_type[1], n_orth = n_orthologs[1])
}, by = symbol]
setorder(spec, panel, -ratio)
## Row label is the marker symbol ONLY (2026-08-20, requested). n_orthologs and the
## specificity ratio still go to marker_specificity.csv -- they are the numbers that decide
## whether a row is usable, they are just no longer printed on the axis.
spec[, row_lab := symbol]
fwrite(spec[order(-ratio)], file.path(OUT, "marker_specificity.csv"))

cat("\n=== marker specificity (best pct / runner-up pct) ===\n")
print(spec[order(-ratio), .(symbol, panel, n_orth, best = round(best, 1),
                            second = round(second, 1), ratio, best_ct)])

d <- merge(d, spec[, .(symbol, row_lab, panel, ratio)], by = "symbol")
d[, panel := factor(panel, levels = unique(spec$panel))]
d[, row_lab := factor(row_lab, levels = rev(unique(spec$row_lab)))]

## Boundaries between tissue groups, in row-index space. row_lab levels are reversed, so
## the break sits between the last row of one panel group and the first of the next.
lv  <- levels(d$row_lab)
grp <- spec[match(lv, row_lab)]$panel
brk <- which(head(grp, -1) != tail(grp, -1)) + 0.5

mk <- function(dd, stem, title, sub, h) {
  n_all <- nrow(dd)
  dd <- dd[expr_scaled >= MIN_SCALED]
  cat(sprintf("%s: kept %d of %d dots at MIN_SCALED = %.2f (%.0f%% hidden)\n",
              stem, nrow(dd), n_all, MIN_SCALED, 100 * (1 - nrow(dd) / n_all)))
  ## Rows whose every dot fell below the cut still get a row -- an empty row means "this
  ## marker is nowhere strong", which is a result, not a gap.
  g <- ggplot(dd, aes(x = col_lab, y = row_lab)) +
    geom_point(aes(size = pct_expressing, fill = expr_scaled),
               shape = 21, colour = "grey40", stroke = 0.25) +
    ## No facets (2026-08-20, requested): 16 strip boxes read as clutter. Row order still
    ## groups markers by the tissue they are reported to mark, and a faint line marks each
    ## boundary, so the grouping survives without the boxes.
    geom_hline(yintercept = brk, colour = "grey75", linewidth = 0.35) +
    scale_fill_gradient(low = LOW, high = HIGH, name = "Scaled mean\nexpression") +
    scale_size_continuous(name = "% nuclei\nexpressing", range = c(0.4, 8),
                          limits = c(0, 100), breaks = c(1, 10, 25, 50, 75)) +
    theme_bw(base_family = FAM, base_size = 12) +
    ## Both axes pinned to the FULL level set. Without this a cluster whose every dot fell
    ## below MIN_SCALED loses its column entirely -- cluster 22 vanished on the first run --
    ## and a reader cannot tell "no strong marker here" from "not in the analysis".
    scale_y_discrete(drop = FALSE, limits = lv) +
    scale_x_discrete(drop = FALSE, limits = levels(d$col_lab)) +
    theme(text = element_text(family = FAM, colour = INK),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = INK, size = 10),
          axis.text.y = element_text(colour = INK, size = 10),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          legend.position = "right") +
    labs(x = NULL, y = NULL, title = title, subtitle = sub)
  f <- file.path(OUT, stem)
  ggsave(paste0(f, ".pdf"), g, width = 13.5, height = h, device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(f, ".png"), g, width = 13.5, height = h, dpi = 300, bg = "white",
         limitsize = FALSE)
  cat(sprintf("wrote %s (%.1f x %.1f in, %d markers)\n", f, 13.5, h, uniqueN(dd$symbol)))
}

SUB <- paste0("Canonical root markers from the literature, grouped by the cell type each is ",
              "reported to mark.\nDot size = % of nuclei expressing; colour = mean ",
              "expression scaled 0-1 within each marker row.\nDots below ",
              MIN_SCALED, " scaled expression are hidden; an empty row means the marker is ",
              "nowhere strong.")

mk(d, "Fig_A_all_canonical_markers",
   "Canonical root markers across res 0.6 clusters", SUB,
   2.8 + 0.185 * uniqueN(d$symbol))

## The publication panel: markers that actually discriminate. A 2x cut is generous -- it
## still admits rows whose best cluster is only twice its runner-up -- but it removes the
## rows that make the full panel look weak (NPF7.3 at 1.2x, ROSY1 1.2x, COBL9 1.5x).
keep <- spec[ratio >= 2]$symbol
cat(sprintf("\nspecific markers (ratio >= 2): %d of %d -- %s\n",
            length(keep), nrow(spec), paste(sort(keep), collapse = ", ")))
mk(d[symbol %in% keep], "Fig_B_specific_markers_only",
   "Canonical markers that discriminate a cluster",
   paste0(SUB, "\nRestricted to markers whose best cluster exceeds its runner-up by >= 2x."),
   2.8 + 0.30 * length(keep))

cat("\ndone:", format(Sys.time()), "\n")
