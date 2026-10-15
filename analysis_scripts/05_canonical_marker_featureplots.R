#!/usr/bin/env Rscript
# Step 05 -- canonical marker FeaturePlots, Control vs Salt side by side,
#            plus the marker table to attach references to.
#
# TWO OUTPUTS
#   1. Fig_canonical_markers_<set>.pdf -- one row per marker, Control | Salt panels
#      side by side on the 500-gate integrated UMAP.
#   2. marker_table_for_reference.csv  -- every marker plotted, with symbol, TAIR
#      locus, the cell type it marks, whether it came from literature or the atlas,
#      the MN106 ortholog actually plotted, how many orthologs it had, and % of
#      nuclei expressing it in Control and in Salt. A blank `reference` column is
#      left for the citation.
#
# ORTHOLOG CHOICE. TAIR symbols map to 1..n MN106 loci (family-level orthology for
# ~77% of genes), so a symbol is NOT a gene here. panel_resolved_orthologs.csv holds
# the OrthoFinder resolution on the same MN106 v4.1 gene space as this object, so it
# is reused rather than recomputed. Where a symbol has several orthologs the one with
# the highest mean expression is plotted, and n_orthologs records how many there
# were -- a marker resolving to 4 loci is weaker evidence than one resolving to 1,
# and the table has to show that rather than hide it.
#
#   MARKER_SET=literature sbatch run_05_marker_featureplots.sh   # classics only (default)
#   MARKER_SET=all        sbatch --export=ALL run_05_marker_featureplots.sh
#
# Outputs -> ncount500_08132026/results/05_marker_featureplots/

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
options(bitmapType = "cairo")

FAM  <- "Times New Roman"
BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
OUT  <- file.path(ROOT, "results/05_marker_featureplots")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PANEL <- file.path(BASE, "results/marker_featureplots_res0.6_073126/marker_panel.csv")
ORTH  <- file.path(BASE, "results/marker_featureplots_res0.6_073126/panel_resolved_orthologs.csv")
stopifnot(file.exists(RDS), file.exists(PANEL), file.exists(ORTH))

MARKER_SET <- Sys.getenv("MARKER_SET", "literature")
REDUC      <- "umap.integrated"

panel <- fread(PANEL)
orth  <- fread(ORTH)
setnames(orth, tolower(names(orth)))
cat("panel markers:", nrow(panel), " ortholog rows:", nrow(orth), "\n")

if (MARKER_SET == "literature") {
  panel <- panel[source == "literature"]
  cat("restricting to literature ('classic') markers:", nrow(panel), "\n")
}

obj <- readRDS(RDS)
stopifnot(REDUC %in% Reductions(obj))
obj$condition <- factor(as.character(obj$condition), levels = c("Control", "Salt"))
DefaultAssay(obj) <- "RNA"
cat("nuclei:", ncol(obj), "\n")

## ---- resolve each symbol to ONE plottable feature ---------------------------
orth <- orth[feature %in% rownames(obj)]
orth <- orth[symbol %in% panel$symbol]
if (!nrow(orth)) stop("no panel ortholog is present in the object -- wrong genome build?")

expr <- GetAssayData(obj, layer = "data")
mean_expr <- Matrix::rowMeans(expr[unique(orth$feature), , drop = FALSE])
orth[, mean_expr := mean_expr[feature]]
setorder(orth, symbol, -mean_expr)

chosen <- orth[, .(feature = feature[1], n_orthologs = .N,
                   mean_expr = round(mean_expr[1], 4)), by = symbol]
chosen <- merge(panel, chosen, by = "symbol", all.x = TRUE)
missing <- chosen[is.na(feature)]
if (nrow(missing)) {
  cat("\nNO ortholog in this object for:", paste(missing$symbol, collapse = ", "), "\n")
  cat("These are dropped from the figure but kept in the table with a blank feature.\n")
}
plotme <- chosen[!is.na(feature)]
setorder(plotme, cell_type, symbol)
cat("\nplotting", nrow(plotme), "markers across", uniqueN(plotme$cell_type), "cell types\n")

## ---- % expressing per condition, for the table ------------------------------
cond <- obj$condition
pct <- rbindlist(lapply(seq_len(nrow(plotme)), function(i) {
  f <- plotme$feature[i]
  v <- expr[f, ] > 0
  data.table(symbol = plotme$symbol[i], feature = f,
             pct_Control = round(100 * mean(v[cond == "Control"]), 1),
             pct_Salt    = round(100 * mean(v[cond == "Salt"]), 1))
}))
tab <- merge(plotme, pct, by = c("symbol", "feature"), all.x = TRUE)
tab[, reference := ""]          # <- fill in the citation
tab <- tab[, .(symbol, tair, cell_type, source, mn106_feature = feature,
               n_orthologs, pct_Control, pct_Salt, reference)]
setorder(tab, cell_type, symbol)
# markers with no ortholog are kept so the table documents the whole panel
if (nrow(missing)) {
  tab <- rbind(tab, missing[, .(symbol, tair, cell_type, source,
                                mn106_feature = NA_character_, n_orthologs = 0L,
                                pct_Control = NA_real_, pct_Salt = NA_real_,
                                reference = "")], fill = TRUE)
}
fwrite(tab, file.path(OUT, "marker_table_for_reference.csv"))
cat("wrote marker table ->", file.path(OUT, "marker_table_for_reference.csv"), "\n")
print(as.data.frame(tab), row.names = FALSE)

## ---- figure: Control | Salt side by side ------------------------------------
# split.by draws one panel per condition on a shared colour scale, which is what
# makes the two directly comparable. keep.scale = "feature" fixes the scale across
# the split so a colour means the same thing in both panels.
th <- theme(text = element_text(family = FAM, colour = "black"),
            axis.text = element_text(colour = "black", size = 6),
            axis.title = element_text(size = 7),
            plot.title = element_text(family = FAM, face = "bold.italic", size = 9),
            legend.text = element_text(size = 6),
            legend.key.width = unit(0.25, "cm"))

make_row <- function(i) {
  f   <- plotme$feature[i]
  ttl <- sprintf("%s  (%s)  %s", plotme$symbol[i], plotme$tair[i],
                 gsub("_", " ", plotme$cell_type[i]))
  FeaturePlot(obj, features = f, reduction = REDUC, split.by = "condition",
              pt.size = 0.15, order = TRUE, combine = TRUE,
              keep.scale = "feature") &
    labs(title = ttl) & th
}

rows <- lapply(seq_len(nrow(plotme)), make_row)
g <- wrap_plots(rows, ncol = 1)

H <- 2.3 * nrow(plotme)
fn <- file.path(OUT, paste0("Fig_canonical_markers_", MARKER_SET))
ggsave(paste0(fn, ".pdf"), g, device = cairo_pdf, width = 8.5, height = H,
       limitsize = FALSE)
ggsave(paste0(fn, ".png"), g, width = 8.5, height = H, dpi = 150, bg = "white",
       limitsize = FALSE)
cat("\nwrote figure ->", paste0(fn, ".pdf"), "\n")

## also a compact per-cell-type version: one page per cell type
for (ctp in unique(plotme$cell_type)) {
  idx <- which(plotme$cell_type == ctp)
  gg  <- wrap_plots(lapply(idx, make_row), ncol = 1)
  ggsave(file.path(OUT, sprintf("Fig_markers_%s.pdf", gsub("[^A-Za-z0-9]+", "_", ctp))),
         gg, device = cairo_pdf, width = 8.5, height = 2.3 * length(idx),
         limitsize = FALSE)
}
cat("wrote per-cell-type figures ->", OUT, "\n")
cat("DONE\n")
