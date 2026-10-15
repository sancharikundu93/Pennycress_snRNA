#!/usr/bin/env Rscript
# The project cell-type palette -- ONE definition, sourced by every figure script in
# this folder, so a cell type is the same colour in the UMAP, the proportion bars and
# anything added later. Copied verbatim from
# troubleshooting_scripts/scripts_07312026/21_umap_split_by_condition.R:56-85, which
# took it from `pastel_cols` in scripts/pericycle_subcluster.Rmd -- the exact colours
# behind results/UMAP/May11_26/UMAP_int_label.pdf.
#
# Pericycle / XPP / PPP / Procambium sit side by side in the embedding, so they take
# four hues ~90 deg apart rather than four shades of one family:
#   Pericycle rose | Procambium violet | XPP dark cyan | PPP dark goldenrod
#
#   source(file.path(dirname(SCRIPTS), "palette_celltypes.R"))
#   -> PAL, CT_LEVELS, shade(), cluster_palette()

PAL <- c(
  "Meristem"              = "#7EC0EE",
  "S phase"               = "#4292C6",
  "G2/M phase"            = "#45B3E0",
  "Root cap"              = "#DA70D6",
  "Trichoblast"           = "#26A69A",
  # lightened 2026-08-13 (was #7F77DD) -- the deep periwinkle read too close to
  # Procambium's violet where the two sit next to each other
  "Atrichoblast"          = "#A49EE8",
  "Epidermis"             = "#6495ED",
  "Cortex"                = "#FFA500",
  # Curated as two distinct cell types 2026-08-13, not as one type split across two
  # clusters: pennycress has three ground tissue layers, so clusters 9 and 14 are
  # taken to be separate cortical layers rather than one population the clustering
  # happened to divide. They therefore carry their own PAL entries and their own
  # rows in every proportion and DE table -- which the automatic _1/_2 suffixing in
  # step 06 never did, since that only renamed a shared label on the plot.
  "Cortex_1"              = "#FFA500",
  "Cortex_2"              = "#C87F00",
  "Endodermis"            = "#41AB5D",
  "Endodermis_1"          = "#77DD77",
  "Endodermis_2"          = "#41AB5D",
  "Endodermis_3"          = "#FFEF00",
  "Pericycle"             = "#D4537E",
  # cluster 5, new at the 500 gate: pericycle rose family but lighter, so it reads as
  # pericycle-adjacent and provisional rather than a fourth confident stele identity
  "pericycle_like"        = "#F0A7C0",
  # lightened 2026-08-13 (was #6A3D9A)
  "Procambium"            = "#9B72CF",
  "Companion cell"        = "#FF85C0",
  "Sieve element"         = "#FFC0CB",
  "Xylem"                 = "#F6B1A7",
  "Xylem pole pericycle"  = "#00838F",
  "Phloem pole pericycle" = "#B8860B",
  "Photosynthetic cell"   = "#9CCC65",
  "Seed/embryo remnant"   = "#D7CCC8",
  # Step 74 re-derived the res 0.6 labels from THIS object with no transfer from the
  # 1000-gate annotation, and renamed two clusters. Both keep their predecessor's colour
  # so every existing figure stays colour-comparable -- only the legend text changes.
  # The old names are KEPT above, not replaced: figures already on disk still carry them.
  "LEA/desiccation-high"     = "#D7CCC8",   # was "Seed/embryo remnant"  (cluster 15)
  "Photosynthetic epidermis" = "#9CCC65")   # was "Photosynthetic cell"  (cluster 22)

# bottom-to-top / outer-to-inner anatomical order, so a stacked bar and a legend both
# read like the root rather than like a frequency table
CT_LEVELS <- c("Meristem", "S phase", "G2/M phase",
               "Root cap", "Trichoblast", "Atrichoblast", "Epidermis",
               "Cortex", "Cortex_1", "Cortex_2", "Endodermis",
               "Phloem pole pericycle", "Pericycle", "pericycle_like",
               "Xylem pole pericycle", "Procambium",
               "Companion cell", "Sieve element", "Xylem",
               "Photosynthetic cell", "Seed/embryo remnant",
               "Photosynthetic epidermis", "LEA/desiccation-high")

## ---- shared UMAP orientation ------------------------------------------------
# rotate_90 chosen from the step 13 comparison, 2026-08-13. UMAP axes carry no
# direction or units, so this is a relabelling of the two stored coordinates:
# distances, neighbours and cluster shapes are untouched. Set here so every figure
# in the folder shares one orientation; override per-run with ORIENT=<mode>.
ORIENT <- Sys.getenv("ORIENT", "rotate_90")

orient_xy <- function(x, y, mode = ORIENT) {
  switch(mode,
         original  = list(x,  y),
         flip_x    = list(-x,  y),
         flip_y    = list(x, -y),
         flip_both = list(-x, -y),
         rotate_90 = list(y, -x),
         stop("unknown ORIENT: ", mode))
}

# Rewrite a reduction in place so every DimPlot downstream inherits the orientation
# without each call having to know about it. The reduction key is preserved, or the
# axis titles change and Seurat starts warning about key mismatches.
orient_reduction <- function(obj, reduc, mode = ORIENT) {
  if (identical(mode, "original")) return(obj)
  e  <- Seurat::Embeddings(obj, reduc)
  xy <- orient_xy(e[, 1], e[, 2], mode)
  e[, 1] <- xy[[1]]; e[, 2] <- xy[[2]]
  obj[[reduc]] <- Seurat::CreateDimReducObject(
    embeddings = e, key = Seurat::Key(obj[[reduc]]),
    assay = Seurat::DefaultAssay(obj))
  obj
}

## ---- on-plot labels ---------------------------------------------------------
# Labels sit over dense points, so they get a white halo behind the glyphs -- the
# "whiteish space" used in scripts_07312026/21_umap_split_by_condition.R:155
# (bg.color = "white", bg.r = 0.16). Size 5 (mm, ~14 pt) and bold, per the
# 2026-08-13 choice.
LABEL_SIZE <- 5

# Seurat names the axes from the reduction key, giving "UMAPintegrated_1". Override
# to plain UMAP_1 / UMAP_2 -- the integration method belongs in the methods text, not
# on every axis of every figure.
AXES <- ggplot2::labs(x = "UMAP_1", y = "UMAP_2")

# cen: a data.frame with x, y and a `label` column (one row per group)
#
# segments = FALSE suppresses the grey leader lines ggrepel draws from a displaced
# label back to its centroid. They read as a stray dash next to the text when the
# displacement is small, which is why the by-replicate figure turns them off; the
# default stays TRUE so 06_labelled_umap is unchanged.
add_labels <- function(p, cen, size = LABEL_SIZE, fam = "Times New Roman",
                       segments = TRUE, bold = TRUE) {
  p + AXES + ggrepel::geom_text_repel(
    data = cen, aes(x = x, y = y, label = label), inherit.aes = FALSE,
    family = fam, size = size, fontface = if (bold) "bold" else "plain",
    colour = "black",
    bg.color = "white", bg.r = 0.16, box.padding = 0.38, point.padding = 0,
    min.segment.length = if (segments) 0.24 else Inf,
    segment.colour = "grey55", segment.size = 0.28,
    max.overlaps = Inf, seed = 42)
}

# lighten (f > 0) or darken (f < 0) a hex colour
shade <- function(hex, f) {
  v <- col2rgb(hex) / 255
  v <- if (f >= 0) v + (1 - v) * f else v * (1 + f)
  rgb(v[1], v[2], v[3])
}

# Colours for a set of labels, erroring on anything the palette does not know rather
# than letting ggplot drop those points to grey -- a silently grey cluster is
# indistinguishable from a real "unassigned".
ct_palette <- function(labels) {
  labels <- as.character(unique(labels))
  miss <- setdiff(labels, names(PAL))
  if (length(miss)) stop("no colour in PAL for: ", paste(miss, collapse = ", "))
  PAL[intersect(CT_LEVELS, labels)]
}

# ---------------------------------------------------------------------------------
# PER-CLUSTER LABELS AND COLOURS: Endodermis_1 / _2 / _3 and their shades.
#
# Lifted out of 06_labelled_umap.R on 2026-08-14, unchanged, because a SECOND figure
# (38_umap_split_by_replicate.R) now needs the same suffixed labels. Step 06's header
# already claimed labels were shared rather than restated; this makes that true. Copying
# the block into the second script instead would have let the two figures disagree about
# which cluster is Endodermis_2, which is exactly the drift the claim warns against.
#
# A cell type spanning several clusters is suffixed _1.._n numbered BY CLUSTER ID -- the
# same convention as scripts_07312026/21_umap_split_by_condition.R and
# results/UMAP/May11_26/UMAP_int_label.pdf. Three identically named regions on one UMAP
# otherwise read as a mistake. k = 1 keeps the canonical unshaded colour; the rest take
# shades either side. Explicit palette entries (Endodermis_1 etc.) win over generated
# shades, which is why the named lookup is tried first.
#
#   cl        character vector of cluster ID per nucleus
#   cell_type character vector of cell type per nucleus, same length
#
# Returns a list: map (cluster/cell_type/label/colour), cl_cols (named by cluster ID,
# ordered numerically), ctn_cols (named by label), ctn_levels, and label_per_cell.
cluster_label_map <- function(cl, cell_type) {
  cl <- as.character(cl); cell_type <- as.character(cell_type)
  stopifnot(length(cl) == length(cell_type))
  present <- sort(unique(cell_type))
  miss <- setdiff(present, names(PAL))
  if (length(miss)) stop("no colour in PAL for: ", paste(miss, collapse = ", "))

  map   <- unique(data.table::data.table(cluster = cl, cell_type = cell_type))
  sizes <- data.table::data.table(cluster = cl)[, .(cells = .N), by = cluster]
  map   <- merge(map, sizes, by = "cluster", all.x = TRUE)
  map[, cluster_num := as.integer(cluster)]
  data.table::setorder(map, cell_type, cluster_num)
  map[, k := seq_len(.N), by = cell_type]
  map[, n := .N, by = cell_type]
  map[, label := ifelse(n > 1, sprintf("%s_%d", cell_type, k), cell_type)]
  map[, colour := {
    base  <- PAL[[cell_type[1]]]
    named <- sprintf("%s_%d", cell_type[1], k)
    fac   <- if (n[1] == 1) 0 else c(0, seq(-0.28, 0.32, length.out = n[1] - 1))
    ifelse(named %in% names(PAL), PAL[named],
           vapply(seq_len(.N), function(i) shade(base, fac[i]), ""))
  }, by = cell_type]

  cl_cols <- setNames(map$colour, map$cluster)
  cl_cols <- cl_cols[levels(factor(cl, levels = sort(as.numeric(unique(cl)))))]
  ctn_levels <- unlist(lapply(intersect(CT_LEVELS, present), function(b) {
    map[cell_type == b][order(cluster_num)]$label
  }), use.names = FALSE)
  # unname(): indexing a named vector by cluster ID returns a vector NAMED by cluster
  # ID, and Seurat's $<- matches names against cell barcodes -- keeping them fails with
  # "No cell overlap between new meta data and Seurat object" (job 668764).
  label_per_cell <- unname(setNames(map$label, map$cluster)[cl])
  stopifnot(!anyNA(label_per_cell))

  list(map = map, cl_cols = cl_cols,
       ctn_cols = setNames(map$colour, map$label)[ctn_levels],
       ctn_levels = ctn_levels, label_per_cell = label_per_cell)
}
