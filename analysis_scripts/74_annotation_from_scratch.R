#!/usr/bin/env Rscript
# =============================================================================
# Step 74 -- re-derive the res 0.6 annotation FROM SCRATCH, with no label transfer
#
#   Rscript 74_annotation_from_scratch.R
#
# WHY (requested 2026-08-19). The labels in 02_annotation/labels_res0.6.csv were not
# derived from this object. 04_annotation_worksheet.py computes three lines of evidence
# and then pre-fills the suggestion from the CROSSTAB PRIOR -- the majority curated label
# of the previous 1000-gate, 18,475-nucleus object (line 70, `crosstab_label`). For
# cluster 15 that prior said "Seed/embryo remnant" at 81.5% support and won, even though:
#
#   * the canonical marker line said Xylem_pole_pericycle (51.1%)
#   * the atlas overlap line said Meristem_Endodermis at FDR 0.148 -- NOT significant
#   * 16.9% of the cluster (82 of 486 nuclei) did not exist in the 1000-gate object at
#     all, so they inherited a label from nuclei they merely clustered with
#
# So the name is a curation decision made once on a different object, three steps removed
# from any measurement on this one. This script drops the crosstab entirely.
#
# THE CANONICAL LINE WAS ALSO BROKEN, AND THAT IS WHY THE PRIOR COULD WIN UNOPPOSED.
# `canonical_best_RNA_snn_res_0_6.csv` scores each panel by MEAN pct_expressing over the
# panel's markers and takes the highest. That ranks panels by how abundantly their genes
# are expressed ANYWHERE, not by how specific they are to a cluster: the
# Xylem_pole_pericycle panel averages 44.7% detection across every cluster in the atlas
# (DOT1 34-59% everywhere, NPF7.3 42-67%), so it wins 16 of 23 clusters and the whole
# line returns only FOUR distinct labels. It is not an annotation, it is a readout of
# panel abundance.
#
# THE FIX, and it is one line of arithmetic: z-score each MARKER across clusters BEFORE
# averaging within a panel. A marker detected at 50% everywhere then contributes nothing;
# a marker at 50% in one cluster and 5% elsewhere dominates. Same input file, same
# numbers, 15 distinct calls instead of 4, and they agree with the curated labels
# wherever the curated label was itself marker-supported (Endodermis c3/c11/c17,
# Root_cap c12, Xylem c19, Root_hair c20, Phloem_companion c8, Procambium c0).
#
# WHAT THIS SCRIPT DOES AND DOES NOT DECIDE.
# It assembles four INDEPENDENT lines of evidence per cluster and writes them side by
# side. It does NOT auto-assign final names, because naming a cluster is a biological
# judgement that needs the marker genes read, and because two of the four lines cannot
# name a cycling cluster at all. Where all lines agree it proposes a call; where they
# disagree it says so and prints the marker genes for curation. Final labels remain a
# human decision -- the point is that the decision is now made against THIS object.
#
# THE FOUR LINES
#   1. canonical_z   specificity-corrected canonical marker panels (the fix above)
#   2. atlas         enrichGO-style overlap against the reference root atlas, with FDR
#                    (02_annotation/annotation_scores_*, unchanged -- this line was fine)
#   3. cellcycle     median S and G2/M module score per CLUSTER. Needed because the two
#                    tissue lines annotate identity, and a cycling cluster's identity is
#                    its cell-cycle phase. Attached by a VALIDATED positional join.
#   4. markers       the top 1-vs-rest markers with TAIR ids, so a human can read them
#
# THE POSITIONAL JOIN, AND WHY IT IS CHECKED.
# 60_cellcycle_module/cellcycle_scores_per_nucleus.csv.gz carries cell_type but NOT
# cluster, and cell_type is the very thing being re-derived -- using it would make the
# line circular. Cluster is recovered from 51_doublets/per_nucleus_doublet_calls.csv.gz
# positionally. Both files are checked to be in the same cell order via their shared
# cell_type column before the join is used, and the script stops if they are not.
#
# Outputs -> ncount500_08132026/results/74_annotation_from_scratch/
# =============================================================================
.libPaths(c("/home/sancharikundu/R/tinkercliffs-rome/4.4.2",
            "/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({ library(data.table) })

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
LIN  <- file.path(BASE, "ncount500_08132026")
RES  <- "0.6"
AN   <- file.path(LIN, "results/02_annotation")
CANP <- file.path(AN, "canonical_marker_pct_RNA_snn_res_0_6.csv")
ATL  <- file.path(AN, "annotation_scores_RNA_snn_res_0_6.csv")
MRK  <- file.path(AN, "markers_TAIR_RNA_snn_res_0_6.csv")
CUR  <- file.path(AN, "labels_res0.6.csv")
CC   <- file.path(LIN, "results/60_cellcycle_module/cellcycle_scores_per_nucleus.csv.gz")
DBL  <- file.path(LIN, "results/51_doublets/per_nucleus_doublet_calls.csv.gz")
OUT  <- file.path(LIN, "results/74_annotation_from_scratch")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("=== step 74: annotation from scratch, res", RES, " (NO label transfer) ===\n")
for (f in c(CANP, ATL, MRK, CUR, CC, DBL)) stopifnot(file.exists(f))

## ---- LINE 1: canonical panels, specificity-corrected -----------------------------------
cp <- fread(CANP)
w  <- dcast(cp, symbol ~ cluster, value.var = "pct_expressing", fun.aggregate = mean)
sym <- w$symbol; mat <- as.matrix(w[, -1]); rownames(mat) <- sym
panel <- unique(cp[, .(symbol, marker_cell_type)])[!duplicated(symbol)]
setkey(panel, symbol)

mu <- rowMeans(mat, na.rm = TRUE)
sd <- apply(mat, 1, sd, na.rm = TRUE)
## A marker with zero variance across clusters carries NO information about which cluster
## it marks. Dropping it is not the same as scoring it 0: leaving it in as 0 would dilute
## its panel's mean toward zero and penalise panels that happen to contain a flat gene.
flat <- which(!is.finite(sd) | sd == 0)
if (length(flat)) cat("dropped", length(flat), "zero-variance canonical markers:",
                      paste(sym[flat], collapse = ", "), "\n")
z <- sweep(sweep(mat, 1, mu, "-"), 1, sd, "/")
if (length(flat)) z <- z[-flat, , drop = FALSE]

pz <- panel[rownames(z), marker_cell_type]
zz <- t(apply(z, 2, function(col) tapply(col, pz, mean, na.rm = TRUE)))
can <- data.table(cluster = as.integer(rownames(zz)),
                  canonical_z      = colnames(zz)[max.col(zz, "first")],
                  canonical_z_score= round(apply(zz, 1, max), 2),
                  canonical_runner = colnames(zz)[apply(zz, 1, function(r) order(-r)[2])],
                  canonical_margin = round(apply(zz, 1, function(r) {
                    s <- sort(r, decreasing = TRUE); s[1] - s[2] }), 2))
fwrite(as.data.table(zz, keep.rownames = "cluster"),
       file.path(OUT, "line1_canonical_z_scores.csv"))
cat("LINE 1: canonical_z gives", uniqueN(can$canonical_z), "distinct calls over",
    nrow(can), "clusters\n")

## ---- LINE 2: reference atlas overlap ---------------------------------------------------
at <- fread(ATL)
atb <- at[order(fdr)][, .SD[1], by = cluster][, .(cluster, atlas = cell_type,
                                                  atlas_fdr = signif(fdr, 3))]
cat("LINE 2: atlas FDR < 0.05 in", sum(atb$atlas_fdr < 0.05), "of", nrow(atb), "clusters\n")

## ---- LINE 3: cell cycle, per CLUSTER, via a validated positional join -------------------
cc  <- fread(CC)
dbl <- fread(DBL)
n1  <- nrow(dbl)
mods <- unique(cc$module)
stopifnot(nrow(cc) == n1 * length(mods))
## CHECK: the two files must share a cell order. Their shared cell_type column is the
## only handle -- if it does not match elementwise for every module block, the join is
## invalid and every cell-cycle number below would be silently wrong.
for (m in mods) {
  blk <- cc[module == m]
  stopifnot(nrow(blk) == n1)
  if (!identical(as.character(blk$cell_type), as.character(dbl$cell_type)))
    stop("CHECK FAILED: cellcycle and doublet files are not in the same cell order ",
         "(module ", m, "). The positional join is invalid.")
}
cat("CHECK passed: cellcycle and doublet files share a cell order (", n1, "nuclei )\n")
cc[, cluster := rep(dbl$cluster, times = length(mods))]
ccw <- dcast(cc[, .(med = median(score)), by = .(cluster, module)],
             cluster ~ module, value.var = "med")
setnames(ccw, setdiff(names(ccw), "cluster"), paste0("cc_", setdiff(names(ccw), "cluster")))
for (j in setdiff(names(ccw), "cluster")) set(ccw, j = j, value = round(ccw[[j]], 3))
fwrite(ccw, file.path(OUT, "line3_cellcycle_by_cluster.csv"))
cat("LINE 3: cell-cycle medians computed for", nrow(ccw), "clusters\n")

## ---- LINE 4: the marker genes, for a human to read -------------------------------------
mk <- unique(fread(MRK), by = c("locus", "cluster"))
top <- mk[order(p_val_adj, -avg_log2FC)][, head(.SD, 12), by = cluster]
tops <- top[, .(top_markers_TAIR = paste(tair, collapse = " ")), by = cluster]
fwrite(top[, .(cluster, locus, tair, avg_log2FC = round(avg_log2FC, 2),
               pct.1, pct.2, p_val_adj)],
       file.path(OUT, "line4_top_markers.csv"))

## ---- assemble --------------------------------------------------------------------------
cur <- fread(CUR)[, .(cluster = as.integer(cluster), current_label = cell_type)]
ev <- Reduce(function(a, b) merge(a, b, by = "cluster", all = TRUE),
             list(cur, can, atb, ccw, tops))

## A cluster is CYCLING if either module median clears the threshold; those are named by
## phase, not tissue, and the two tissue lines cannot speak to them.
CC_CUT <- 0.10
ccols <- grep("^cc_", names(ev), value = TRUE)
ev[, cycling := Reduce(`|`, lapply(ccols, function(cc_) ev[[cc_]] > CC_CUT))]

## Agreement is deliberately CONSERVATIVE. It requires the atlas call to be significant
## AND the canonical call to be a substring-compatible match, because the two vocabularies
## differ (atlas "Maturation_Endodermis" vs canonical "Endodermis"). Anything else is
## flagged for curation rather than silently resolved -- the whole failure being corrected
## here came from silently resolving a disagreement.
norm <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
ev[, agree := mapply(function(a, c_) {
  if (is.na(a) || is.na(c_)) return(FALSE)
  grepl(norm(c_), norm(a), fixed = TRUE) || grepl(norm(a), norm(c_), fixed = TRUE)
}, atlas, canonical_z)]
ev[, atlas_sig := atlas_fdr < 0.05]
ev[, proposed := fifelse(cycling, "CYCLING - name by phase",
                  fifelse(agree & atlas_sig, canonical_z, "CURATE - lines disagree"))]

setcolorder(ev, c("cluster", "current_label", "proposed", "canonical_z",
                  "canonical_z_score", "canonical_margin", "canonical_runner",
                  "atlas", "atlas_fdr", ccols, "cycling", "agree"))
fwrite(ev, file.path(OUT, "annotation_evidence_res0.6.csv"))

cat("\n=== EVIDENCE, no label transfer anywhere ===\n")
print(ev[, .(cluster, current_label, canonical_z, canonical_margin, atlas, atlas_fdr,
             cycling, proposed)])

cat("\n=== clusters whose current label has NO support from either tissue line ===\n")
bad <- ev[!cycling & !agree]
for (i in seq_len(nrow(bad))) {
  r <- bad[i]
  cat(sprintf("\ncluster %-2s  current: %-22s\n  canonical_z: %-22s (z %.2f, margin %.2f)\n  atlas:       %-22s (FDR %s)\n  markers:     %s\n",
              r$cluster, r$current_label, r$canonical_z, r$canonical_z_score,
              r$canonical_margin, r$atlas, format(r$atlas_fdr), r$top_markers_TAIR))
}

cat("\nwrote ->", OUT, "\n")
cat("done:", format(Sys.time()), "\n")
