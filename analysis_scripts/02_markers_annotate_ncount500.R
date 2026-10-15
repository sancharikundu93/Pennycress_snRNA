#!/usr/bin/env Rscript
# Step 02 -- markers per cluster + annotation EVIDENCE for the 500-gate object.
#
# This produces the evidence for annotation. It does NOT assign final cell-type
# names -- review the outputs and set the labels yourself. The earlier annotation's
# ~30% pericycle with zero procambium is exactly the kind of call that should not be
# made automatically.
#
# Mirrors step_06/02_markers_annotate.R (same markers, same TAIR mapping, same
# hypergeometric against the root atlas) and adds ONE thing that step_06 could not
# have: a cross-tab of every new cluster against the curated labels of the 1000-gate
# object. 18,475 of these 25,206 nuclei are in both objects and carry a curated
# label there, so most of the 23 new clusters should resolve to an obvious majority
# label rather than being annotated from scratch. Treat it as a strong prior, not an
# answer -- the 3 extra clusters and any split/merged ones need real marker review.
#
#   CLUSTER_COL=RNA_snn_res.0.6 sbatch run_02_markers_annotate.sh
#
# Outputs -> ncount500_08132026/results/02_annotation/

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(data.table)
})

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
OUT  <- file.path(ROOT, "results/02_annotation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# The alternative_splicing tree was moved off songli_lab to the intro2gds fileset
# (quota), so the path baked into step_06/02_markers_annotate.R is dead -- it points
# at BASE/alternative_splicing/... which no longer exists. Job 668198 died here.
MAP   <- "/projects/intro2gds/A_Sanchari/alternative_splicing/GO/mapping/tamn106_to_tair.csv"
ATLAS <- file.path(BASE, "markergene/ATh_marker_atlas_rachel_sk.csv")
stopifnot(file.exists(MAP), file.exists(ATLAS))

CLUSTER_COL <- Sys.getenv("CLUSTER_COL", "RNA_snn_res.0.6")
tag <- gsub("[^A-Za-z0-9]+", "_", CLUSTER_COL)

# The cross-tab reference must MATCH THE RESOLUTION BEING ANNOTATED. This used to be
# hardcoded to the old res 0.6 labels, which is right only when CLUSTER_COL is res
# 0.6 -- crosstabbing res 0.25 clusters against a res 0.6 annotation would import
# that resolution's splitting decisions into a coarser one and make the suggestion
# column quietly wrong. marker_panel_root carries curated labels for all five
# resolutions, so use the one that corresponds.
res_num <- sub("^RNA_snn_res\\.", "", CLUSTER_COL)
old_res <- if (identical(res_num, "1")) "1.0" else res_num
OLD_RDS <- file.path(BASE, "results/cluster_proportion_mn106/penny_mn106_integrated.rds")
OLD_LAB <- file.path(BASE, "results/marker_panel_root",
                     sprintf("labels_res%s.csv", old_res))
OLD_COL <- CLUSTER_COL
if (!file.exists(OLD_LAB))
  stop("no curated 1000-gate labels for resolution ", old_res, ": ", OLD_LAB)
cat("cross-tab reference:", basename(OLD_LAB), "  column:", OLD_COL, "\n")

# canonical marker panel -- the third independent line of evidence, alongside the
# atlas hypergeometric and the cross-tab
PANEL <- file.path(BASE, "results/marker_featureplots_res0.6_073126/marker_panel.csv")
ORTH  <- file.path(BASE, "results/marker_featureplots_res0.6_073126/panel_resolved_orthologs.csv")

cat("cluster column:", CLUSTER_COL, "\n")
so <- readRDS(RDS)
stopifnot(CLUSTER_COL %in% colnames(so@meta.data))
Idents(so) <- CLUSTER_COL
cat("nuclei:", ncol(so), " clusters:", length(unique(Idents(so))), "\n")

cat("\n########## 1. markers per cluster ##########\n")
# cached: FindAllMarkers is slow and does not change when only the scoring below is
# edited. Delete markers_<tag>.csv to force recomputation.
mk_file <- file.path(OUT, paste0("markers_", tag, ".csv"))
if (file.exists(mk_file)) {
  mk <- as.data.frame(fread(mk_file))
  cat("reusing cached markers:", mk_file, "\n")
} else {
  mk <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.5, test.use = "wilcox")
  fwrite(mk, mk_file)
}
cat("markers:", nrow(mk), "across", length(unique(mk$cluster)), "clusters\n")

cat("\n########## 2. map to Arabidopsis ##########\n")
# tamn106_to_tair.csv holds bare loci; our features carry the .v4.1 suffix
map <- fread(MAP)
setnames(map, c("mn106_locus", "tair"))
map[, tair := toupper(tair)]
mk$locus <- sub("\\.v4\\.1$", "", mk$gene)
mkm <- merge(as.data.table(mk), map, by.x = "locus", by.y = "mn106_locus",
             allow.cartesian = TRUE)
cat(sprintf("markers mapped to TAIR: %d of %d (%.0f%% of loci)\n",
            nrow(mkm), nrow(mk), 100 * mean(unique(mk$locus) %in% map$mn106_locus)))
fwrite(mkm, file.path(OUT, paste0("markers_TAIR_", tag, ".csv")))

cat("\n########## 3. score against the Arabidopsis root atlas ##########\n")
atlas <- fread(ATLAS)
stopifnot(all(c("cluster", "gene") %in% names(atlas)))
atlas[, gene := toupper(gene)]
# split FIRST, then unique WITHIN each group. split(unique(x), f) silently recycles
# because unique(x) is shorter than f, which scrambles the gene sets.
atlas_sets <- lapply(split(atlas$gene, atlas$cluster), unique)
cat("atlas cell types:", length(atlas_sets), "\n")

# background = every TAIR gene reachable from our expressed feature space, so the
# hypergeometric is conditioned on what could have been detected at all
bg <- unique(map[mn106_locus %in% sub("\\.v4\\.1$", "", rownames(so)), tair])
cat("background TAIR genes:", length(bg), "\n\n")

ours <- lapply(split(mkm$tair, mkm$cluster), unique)
res <- rbindlist(lapply(names(ours), function(cl) {
  q <- intersect(ours[[cl]], bg)
  rbindlist(lapply(names(atlas_sets), function(ct) {
    s <- intersect(atlas_sets[[ct]], bg)
    ov <- length(intersect(q, s))
    p <- if (ov == 0) 1 else phyper(ov - 1, length(s), length(bg) - length(s),
                                    length(q), lower.tail = FALSE)
    data.table(cluster = cl, cell_type = ct, n_markers = length(q),
               n_atlas = length(s), overlap = ov, p = p)
  }))
}))
res[, fdr := p.adjust(p, method = "BH")]
setorder(res, cluster, p)
fwrite(res, file.path(OUT, paste0("annotation_scores_", tag, ".csv")))

cat("top 3 atlas matches per cluster:\n")
print(as.data.frame(res[, head(.SD, 3), by = cluster][
  , .(cluster, cell_type, overlap, n_markers, fdr)]), digits = 3)

cat("\n########## 4. cross-tab against the 1000-gate curated labels ##########\n")
new_meta <- data.table(cell = colnames(so),
                       new_cluster = as.character(so@meta.data[[CLUSTER_COL]]))
sz <- new_meta[, .(cells = .N), by = new_cluster]
sz[, pct := round(100 * cells / sum(cells), 1)]
## ---- canonical marker panel: % of nuclei expressing, per cluster -------------
# The third independent line of evidence. The atlas hypergeometric asks "is this
# cluster's marker LIST enriched for a reference cell type"; this asks the blunter
# question "do the canonical markers of that cell type actually fire in these
# nuclei". They fail in different ways, which is why both are kept.
#
# A TAIR symbol maps to 1..n MN106 loci (family level for most genes), so the
# highest-expressed ortholog is used and n_orthologs is carried through -- a 9-
# ortholog "marker" is not gene-level evidence and the table must say so.
if (file.exists(PANEL) && file.exists(ORTH)) {
  cat("\n########## canonical marker panel ##########\n")
  panel <- fread(PANEL)
  orth  <- fread(ORTH); setnames(orth, tolower(names(orth)))
  orth  <- orth[feature %in% rownames(so)][symbol %in% panel$symbol]
  if (nrow(orth)) {
    expr <- GetAssayData(so, layer = "data")
    me   <- Matrix::rowMeans(expr[unique(orth$feature), , drop = FALSE])
    orth[, mean_expr := me[feature]]
    setorder(orth, symbol, -mean_expr)
    chosen <- orth[, .(feature = feature[1], n_orthologs = .N), by = symbol]
    chosen <- merge(panel, chosen, by = "symbol")[!is.na(feature)]

    cl_vec <- as.character(so@meta.data[[CLUSTER_COL]])
    m <- expr[chosen$feature, , drop = FALSE]
    pct <- do.call(rbind, lapply(sort(unique(cl_vec)), function(k) {
      i <- which(cl_vec == k)
      data.table(cluster = k, symbol = chosen$symbol,
                 marker_cell_type = chosen$cell_type,
                 n_orthologs = chosen$n_orthologs,
                 pct_expressing = round(100 * Matrix::rowSums(m[, i, drop = FALSE] > 0) /
                                        length(i), 2),
                 mean_expr = round(Matrix::rowMeans(m[, i, drop = FALSE]), 4))
    }))
    fwrite(pct, file.path(OUT, paste0("canonical_marker_pct_", tag, ".csv")))

    # per cluster, which marker cell type has the highest mean detection
    top <- pct[, .(mean_pct = mean(pct_expressing)), by = .(cluster, marker_cell_type)]
    setorder(top, cluster, -mean_pct)
    best <- top[, .(canonical_best = marker_cell_type[1],
                    canonical_pct = round(mean_pct[1], 2),
                    runner_up = marker_cell_type[2],
                    runner_up_pct = round(mean_pct[2], 2)), by = cluster]
    fwrite(best, file.path(OUT, paste0("canonical_best_", tag, ".csv")))
    cat("canonical panel: ", nrow(chosen), " markers over ", uniqueN(pct$cluster),
        " clusters\n", sep = "")
    print(best)
  } else cat("no panel orthologs present in this object -- skipped\n")
} else cat("\ncanonical panel files absent -- skipped\n")

rm(so); invisible(gc())

if (file.exists(OLD_RDS) && file.exists(OLD_LAB)) {
  old <- readRDS(OLD_RDS)
  old_meta <- data.table(cell = colnames(old),
                         old_cluster = as.character(old@meta.data[[OLD_COL]]))
  rm(old); invisible(gc())
  lab <- fread(OLD_LAB, select = c("cluster", "cell_type"))
  lab[, cluster := as.character(cluster)]
  old_meta[, old_label := lab$cell_type[match(old_cluster, lab$cluster)]]

  j <- merge(new_meta, old_meta[, .(cell, old_label)], by = "cell", all.x = TRUE)
  j[is.na(old_label), old_label := "(not in 1000-gate object)"]
  cat("nuclei shared with the 1000-gate object:",
      sum(j$old_label != "(not in 1000-gate object)"), "of", nrow(j), "\n\n")

  ct <- j[, .N, by = .(new_cluster, old_label)]
  ct[, pct_of_new := round(100 * N / sum(N), 1), by = new_cluster]
  setorder(ct, new_cluster, -N)
  fwrite(ct, file.path(OUT, paste0("crosstab_new_vs_old_", tag, ".csv")))

  # majority old label per new cluster, ignoring nuclei the 1000 gate dropped
  top <- ct[old_label != "(not in 1000-gate object)"][
    order(new_cluster, -N)][, head(.SD, 1), by = new_cluster]
  sug <- merge(sz, top[, .(new_cluster, suggested = old_label,
                           support_pct = pct_of_new)], by = "new_cluster", all.x = TRUE)
  best <- res[order(cluster, fdr)][, head(.SD, 1), by = cluster][
    , .(new_cluster = cluster, atlas_best = cell_type, atlas_fdr = signif(fdr, 3))]
  sug <- merge(sug, best, by = "new_cluster", all.x = TRUE)
  sug[, new_cluster := factor(new_cluster, levels = as.character(sort(as.numeric(new_cluster))))]
  setorder(sug, new_cluster)
  fwrite(sug, file.path(OUT, paste0("label_suggestions_", tag, ".csv")))

  cat("SUGGESTED labels -- majority 1000-gate label vs best atlas match.\n")
  cat("Where the two columns disagree, review the markers before deciding.\n\n")
  print(as.data.frame(sug), digits = 3)
} else {
  cat("SKIPPED -- 1000-gate object or labels not found:\n  ", OLD_RDS, "\n  ", OLD_LAB, "\n")
  fwrite(sz, file.path(OUT, paste0("cluster_sizes_", tag, ".csv")))
}

cat("\nNEXT: curate the labels into labels_", tag, ".csv (cluster,cell_type),\n", sep = "")
cat("then run 03_proportions.\n")
cat("DONE\n")
