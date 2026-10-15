#!/usr/bin/env Rscript
# =============================================================================
# Step 77 -- the annotation-justification dot plot: TWO specific markers per cluster
#
#   sbatch run_77_top2_dotplot.sh
#
# WHY NOT THE CANONICAL PANEL (step 76). The literature panel has 44 markers and only TEN
# clear a 2x specificity bar, and those ten do not cover the atlas -- they cluster on
# endodermis, trichoblast and xylem and leave most clusters with nothing. A figure that
# justifies an annotation has to speak for EVERY cluster, so the markers are taken from
# this object's own 1-vs-rest FindMarkers table instead of from a fixed literature list.
# The canonical panel remains the independent check (step 76); this is the descriptive one.
#
# HOW THE TWO PER CLUSTER ARE CHOSEN. Candidates are filtered to pct.1 > 0.30,
# avg_log2FC > 1, p_val_adj < 1e-10, then ranked by the DETECTION GAP pct.1 - pct.2.
# The gap, not the fold change: a gene can post a huge log2FC off a handful of counts, but
# a large gap means the gene is on in most of this cluster and off in most of the rest,
# which is exactly the claim a marker figure makes. Each locus is then assigned to the ONE
# cluster it separates best, so no gene is claimed twice and the diagonal is real rather
# than an artefact of showing each gene against its own cluster.
#
# DISTINCT TAIR PER CLUSTER is enforced. Without it cluster 22 draws TaMN106.1G095800 and
# .1G095900 -- two loci that both best-hit AT1G29910 (CAB/LHCB1) -- and the row pair looks
# like two independent markers when it is one Arabidopsis gene counted twice.
#
# WHY THIS ONE NEEDS THE OBJECT AND STEP 76 DID NOT. A dot plot needs each gene's detection
# in EVERY cluster. canonical_marker_pct_*.csv already stored that for the literature panel;
# for arbitrary genes it has to be computed from the counts matrix. Job 696963 showed the
# 5 GB object loads in about a minute, so this is cheap despite the memory request.
#
# WHAT TO READ OFF IT. Dot size = % of nuclei expressing; colour = mean expression scaled
# 0-1 within each row. A convincing panel shows one large dark dot per row sitting on the
# diagonal, and small pale dots everywhere else. Rows that fail that test are labelled with
# their specificity ratio so the failure is visible rather than argued about.
#
# Outputs -> ncount500_08132026/results/77_top2_dotplot/
# =============================================================================
.libPaths(c("/home/sancharikundu/R/tinkercliffs-rome/4.4.2",
            "/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); library(Matrix)
  library(org.At.tair.db); library(AnnotationDbi)
})
options(bitmapType = "cairo")

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
LIN  <- file.path(BASE, "ncount500_08132026")
RES  <- "0.6"
RDS  <- file.path(LIN, "results/penny_mn106_ncount500_integrated.rds")
MRK  <- file.path(LIN, "results/02_annotation/markers_TAIR_RNA_snn_res_0_6.csv")
## Labels reverted to the ORIGINAL file (2026-08-20, requested): "Photosynthetic cell" and
## "Seed/embryo remnant" stay for now. Step 74's from-scratch table still exists and still
## disagrees with both names -- this is a presentation choice, not a retraction of that
## analysis. Point LABELS at the step-74 file to switch back.
LAB  <- Sys.getenv("LABELS", file.path(LIN, "results/02_annotation/labels_res0.6.csv"))
OUT  <- file.path(LIN, "results/77_top2_dotplot")
source(file.path(LIN, "scripts/palette_celltypes.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FAM <- "Times New Roman"; INK <- "black"
LOW <- "#F2F0F7"; HIGH <- "#4A1486"
N_PER <- as.integer(Sys.getenv("N_PER", "2"))
## Clusters dropped from the figure entirely -- rows AND columns (2026-08-20, requested).
## 5 = pericycle_like and 15 = Seed/embryo remnant are the two whose identity is still
## open (step 74: both flagged, specificity 1.1-1.2x and the marker set disputed), so they
## are held out of the annotation-justification panel rather than shown as settled.
DROP <- as.integer(strsplit(Sys.getenv("DROP_CLUSTERS", "5,15"), ",")[[1]])
cat("dropping clusters:", paste(DROP, collapse = ", "), "\n")

cat("=== step 77: top", N_PER, "markers per cluster, res", RES, "===\n")
stopifnot(file.exists(RDS), file.exists(MRK), file.exists(LAB))

## ---- choose the markers ----------------------------------------------------------------
lab <- fread(LAB)[, .(cluster = as.integer(cluster), cell_type)]
lab <- lab[!cluster %in% DROP]
m   <- unique(fread(MRK), by = c("locus", "cluster"))
m    <- m[!cluster %in% DROP]
cand <- m[`pct.1` > 0.30 & avg_log2FC > 1 & p_val_adj < 1e-10]
cand[, dpct := `pct.1` - `pct.2`]
setorder(cand, -dpct)
cand <- unique(cand, by = "locus")             # each locus -> the cluster it separates best
## distinct Arabidopsis gene per cluster, so a tandem pair does not fill both slots
cand <- unique(cand, by = c("cluster", "tair"))
sel <- cand[, head(.SD, N_PER), by = cluster]
sel <- merge(sel, lab, by = "cluster")
cat("selected", nrow(sel), "markers over", uniqueN(sel$cluster), "clusters\n")
if (uniqueN(sel$cluster) < uniqueN(lab$cluster))
  warning("clusters with no marker passing the filters: ",
          paste(setdiff(lab$cluster, sel$cluster), collapse = ", "))

## ---- object: detection and mean expression per cluster ---------------------------------
cat("loading object ...\n")
so <- readRDS(RDS)
res_col <- paste0("RNA_snn_res.", RES)
stopifnot(res_col %in% colnames(so@meta.data))
cl <- as.integer(as.character(so@meta.data[[res_col]]))
cat("object:", ncol(so), "nuclei x", nrow(so), "genes\n")

feat <- paste0(sel$locus, ".v4.1")
miss <- setdiff(feat, rownames(so))
if (length(miss)) stop("not in the matrix: ", paste(miss, collapse = ", "))
X <- GetAssayData(so, assay = "RNA", layer = "data")[feat, , drop = FALSE]

cls <- sort(unique(cl))
pct  <- sapply(cls, function(k) Matrix::rowMeans(X[, cl == k, drop = FALSE] > 0) * 100)
mexp <- sapply(cls, function(k) Matrix::rowMeans(X[, cl == k, drop = FALSE]))
colnames(pct) <- colnames(mexp) <- cls
rm(so, X); invisible(gc())

d <- merge(
  melt(as.data.table(pct, keep.rownames = "feature"), id.vars = "feature",
       variable.name = "cluster", value.name = "pct"),
  melt(as.data.table(mexp, keep.rownames = "feature"), id.vars = "feature",
       variable.name = "cluster", value.name = "mean_expr"),
  by = c("feature", "cluster"))
d[, cluster := as.integer(as.character(cluster))]
sel[, feature := paste0(locus, ".v4.1")]
d <- merge(d, sel[, .(feature, own_cluster = cluster, tair, locus)], by = "feature",
           allow.cartesian = TRUE)
d <- merge(d, lab, by = "cluster")

## colour scaled 0-1 within a row, so rows are comparable regardless of absolute abundance
d[, expr_scaled := { r <- range(mean_expr); if (diff(r) == 0) 0 else (mean_expr - r[1]) / diff(r) },
  by = feature]

## ---- ordering, and the specificity number that goes on each row label -------------------
ct_ord <- CT_LEVELS[CT_LEVELS %in% unique(lab$cell_type)]
lab[, ct_rank := match(cell_type, ct_ord)]
setorder(lab, ct_rank, cluster)
lab[, col_lab := sprintf("%d  %s", cluster, cell_type)]
d <- merge(d, lab[, .(cluster, col_lab)], by = "cluster")
d[, col_lab := factor(col_lab, levels = lab$col_lab)]

sp <- d[, { s <- sort(pct, decreasing = TRUE)
            .(own = pct[cluster == own_cluster], ratio = round(s[1] / max(s[2], 0.1), 1),
              rank_own = which(order(-pct) == which(cluster == own_cluster))) },
        by = .(feature, tair, own_cluster)]
sp <- merge(sp, lab[, .(own_cluster = cluster, ct_rank)], by = "own_cluster")
setorder(sp, ct_rank, -own)
## ---- AGI -> gene symbol -----------------------------------------------------------------
## Symbols come from org.At.tair.db, not from memory. Most AGIs carry several aliases
## (AT5G57090 is AGR/AGR1/ATPIN2/EIR1/MM31/PIN2/WAV6) and no rule picks the conventional
## one reliably, so the handful that appear in this panel are pinned by name and everything
## else falls back to: drop AT-prefixed aliases, then take the shortest. An AGI with no
## symbol at all keeps its accession -- inventing a name would be worse than showing one.
PIN <- c(AT5G57090 = "PIN2",   AT1G32450 = "NRT1.5", AT3G23430 = "PHO1",
         AT1G02950 = "GSTF4",  AT1G66200 = "GLN1;2", AT3G15353 = "MT3",
         AT4G19840 = "PP2-A1", AT1G29910 = "LHCB1.2", AT1G52400 = "BGLU18",
         AT3G20470 = "GRP5",   AT2G05440 = "GRP9",   AT1G28330 = "DRM1",
         AT5G13870 = "XTH5",   AT5G09530 = "PELPK1", AT5G41040 = "ASFT",
         AT5G24800 = "bZIP9",  AT1G75500 = "WAT1",   AT3G17520 = "LEA30",
         AT5G59870 = "HTA6",   AT3G53740 = "RPL36",  AT1G67090 = "RBCS1A")
agi <- unique(sp$tair)
ann <- tryCatch(as.data.table(AnnotationDbi::select(org.At.tair.db, keys = agi,
                  keytype = "TAIR", columns = "SYMBOL")),
                error = function(e) data.table(TAIR = character(), SYMBOL = character()))
ann <- ann[!is.na(SYMBOL), .(alias = list(unique(SYMBOL))), by = TAIR]
pick <- function(a) {
  a <- a[nzchar(a)]
  if (!length(a)) return(NA_character_)
  b <- a[!grepl("^AT[A-Z]", a)]              # drop ATxxx-style aliases when possible
  if (length(b)) a <- b
  a[which.min(nchar(a))]
}
ann[, sym_auto := vapply(alias, pick, character(1))]
sp <- merge(sp, ann[, .(tair = TAIR, sym_auto)], by = "tair", all.x = TRUE)
sp[, gene := fifelse(tair %in% names(PIN), PIN[tair],
              fifelse(!is.na(sym_auto), sym_auto, tair))]
nosym <- sp[gene == tair]
if (nrow(nosym)) cat("no symbol available, showing the accession:",
                     paste(sort(unique(nosym$tair)), collapse = ", "), "\n")
setorder(sp, ct_rank, -own)
## Row label is the gene name ONLY (2026-08-20, requested): no cluster tag, no ratio.
## Both are still written to selected_marker_specificity.csv.
sp[, row_lab := gene]
sp[, row_lab := make.unique(row_lab, sep = " ")]   # two loci can share one symbol
fwrite(sp[order(ct_rank, -own)], file.path(OUT, "selected_marker_specificity.csv"))
cat("\n=== selected markers, specificity = own cluster %% / next best %% ===\n")
print(sp[order(ct_rank, -own), .(gene, tair, own_cluster, own = round(own, 1), ratio)])

d <- merge(d, sp[, .(feature, row_lab)], by = "feature")
d[, row_lab := factor(row_lab, levels = rev(unique(sp$row_lab)))]
fwrite(d[order(row_lab, col_lab)], file.path(OUT, "dotplot_values_res0.6.csv"))

g <- ggplot(d, aes(x = col_lab, y = row_lab)) +
  geom_point(aes(size = pct, fill = expr_scaled), shape = 21, colour = "grey40",
             stroke = 0.25) +
  scale_fill_gradient(low = LOW, high = HIGH, name = "Scaled mean\nexpression") +
  scale_size_continuous(name = "% nuclei\nexpressing", range = c(0.3, 7.5),
                        limits = c(0, 100), breaks = c(1, 10, 25, 50, 75)) +
  theme_bw(base_family = FAM, base_size = 12) +
  theme(text = element_text(family = FAM, colour = INK),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = INK, size = 10),
        axis.text.y = element_text(colour = INK, size = 9.5),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10)) +
  labs(x = NULL, y = NULL,
       title = sprintf("Two cluster-specific markers per cluster (resolution %s)", RES),
       subtitle = paste0("Dot size = % of nuclei expressing; colour = mean expression ",
                         "scaled 0-1 within each row.\nRow label is the Arabidopsis ",
                         "symbol of the best hit (accession shown where TAIR has no ",
                         "symbol).\nMarkers are the two with the largest detection gap ",
                         "per cluster, from this object's own 1-vs-rest test."))

W <- 5.0 + 0.42 * uniqueN(d$col_lab)
H <- 3.0 + 0.26 * uniqueN(d$row_lab)
f <- file.path(OUT, sprintf("Fig_top%d_marker_dotplot_res%s", N_PER, RES))
ggsave(paste0(f, ".pdf"), g, width = W, height = H, device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(f, ".png"), g, width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("\nwrote %s (%.1f x %.1f in)\n", f, W, H))
cat("done:", format(Sys.time()), "\n")
