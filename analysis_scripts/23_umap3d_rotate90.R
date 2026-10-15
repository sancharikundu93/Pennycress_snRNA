#!/usr/bin/env Rscript
# Step 23 -- 3D UMAP, oriented to match the rotate_90 2D figure.
#
# READ THIS BEFORE COMPARING IT TO THE 2D FIGURE. A 3D UMAP is NOT the 2D one seen
# from a different angle. UMAP optimises the embedding for the number of dimensions
# it is given, so n.components = 3 is a SEPARATE run: cluster positions, and which
# clusters end up adjacent, will differ from umap.integrated. Nothing upstream
# changes -- same integrated.cca reduction, same 30 dims, same seed -- but the
# embedding is new and the 2D figure remains the one the paper's claims rest on.
#
# The rotate_90 convention (palette_celltypes.R, chosen from the step 13 comparison)
# is a relabelling of two stored coordinates, so it carries over to the first two
# axes here; the third is left as UMAP returns it. That keeps the front view reading
# in the same direction as the published 2D panel.
#
# WHAT IT WRITES
#   1. umap3d_embedding.csv  -- the coordinates, per nucleus, with cell type and
#      condition. Small, and it means nobody has to re-run UMAP to replot.
#   2. Fig_umap3d_angles.pdf/.png -- a grid of static views. The projection is done
#      here with an explicit rotation matrix rather than by a 3D library, so the
#      panels are ordinary ggplot output: vector PDF, same fonts, same palette,
#      no rgl/X11 dependency on a compute node.
#   3. Fig_umap3d_<view>.png -- each view on its own, full size.
#   4. umap3d_interactive.html -- rotatable, IF plotly and htmlwidgets both load.
#      htmlwidgets is absent from the current R_LIBS_USER, so this is expected to
#      skip; the static panels are the deliverable and the HTML is a convenience.
#
# Depth cueing: points are drawn back-to-front and far points are smaller and more
# transparent. Without it an orthographic 3D scatter reads as a flat smear.
#
#   sbatch run_23_umap3d.sh
#   VIEWS=8 sbatch --export=ALL run_23_umap3d.sh
#
# Outputs -> ncount500_08132026/results/23_umap3d/
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
options(bitmapType = "cairo")

FAM  <- "Times New Roman"
BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
LAB  <- file.path(ROOT, "results/02_annotation/labels_res0.6.csv")
OUT  <- file.path(ROOT, "results/23_umap3d")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/palette_celltypes.R"))
stopifnot(file.exists(RDS), file.exists(LAB))

CLUSTER_COL <- Sys.getenv("CLUSTER_COL", "RNA_snn_res.0.6")
DIMS        <- 1:30
SEED        <- 42

obj <- readRDS(RDS)
lab <- fread(LAB, colClasses = "character")
cl  <- as.character(obj@meta.data[[CLUSTER_COL]])
hit <- match(cl, lab$cluster)
if (anyNA(hit)) stop("unlabelled cluster(s): ",
                     paste(sort(unique(cl[is.na(hit)])), collapse = ", "))
obj$cell_type <- factor(lab$cell_type[hit],
                        levels = CT_LEVELS[CT_LEVELS %in% lab$cell_type])
cat("nuclei:", ncol(obj), " cell types:", nlevels(obj$cell_type), "\n")

## ---- the 3D embedding -------------------------------------------------------
# Same input reduction and dims as the 2D run in 01_build_cluster_ncount500.R:127-129,
# so the only difference is n.components.
obj <- RunUMAP(obj, dims = DIMS, reduction = "integrated.cca", n.components = 3,
               seed.use = SEED, reduction.name = "umap3d",
               reduction.key = "UMAP3d_", verbose = FALSE)
e <- Embeddings(obj, "umap3d")
cat("embedding:", nrow(e), "x", ncol(e), "\n")

# rotate_90 on the first two axes, exactly as orient_xy() defines it, so the front
# view matches the published 2D orientation. The third axis is untouched.
xy <- orient_xy(e[, 1], e[, 2], "rotate_90")
d <- data.table(x = xy[[1]], y = xy[[2]], z = e[, 3],
                cell_type = obj$cell_type,
                condition = obj$condition,
                cluster   = cl)
fwrite(d, file.path(OUT, "umap3d_embedding.csv"))

## ---- orthographic projection ------------------------------------------------
# Standard camera: azimuth a spins the scene about the vertical axis, elevation e
# tilts it. depth is the along-view distance, used only for draw order and cueing.
project <- function(d, a_deg, e_deg) {
  a <- a_deg * pi / 180; el <- e_deg * pi / 180
  h <- d$x * cos(a) + d$y * sin(a)          # in-view horizontal component
  data.table(u     = -d$x * sin(a) + d$y * cos(a),
             v     = -h * sin(el) + d$z * cos(el),
             depth =  h * cos(el) + d$z * sin(el),
             cell_type = d$cell_type)
}

PAL_CT <- ct_palette(levels(d$cell_type))

mk_view <- function(a_deg, e_deg, title, legend = FALSE) {
  p <- project(d, a_deg, e_deg)
  setorder(p, depth)                         # far points first, near points on top
  rng <- range(p$depth)
  p[, w := if (diff(rng) > 0) (depth - rng[1]) / diff(rng) else 0.5]
  g <- ggplot(p, aes(u, v, colour = cell_type)) +
    # First pass used size 0.10-0.32 and alpha 0.35-0.80: the far half of the cloud
    # washed out to nothing and the panel read as empty. The cue only has to be
    # enough to tell front from back, so the floor is now high enough that a distant
    # nucleus is still a visible point.
    geom_point(aes(size = 0.30 + 0.30 * w, alpha = 0.60 + 0.35 * w),
               stroke = 0, shape = 16) +
    scale_colour_manual(values = PAL_CT, name = "Cell type",
                        guide = if (legend) guide_legend(override.aes =
                          list(size = 3, alpha = 1)) else "none") +
    scale_size_identity() + scale_alpha_identity() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_family = FAM) +
    theme(plot.title = element_text(family = FAM, size = 10, hjust = 0.5,
                                    colour = "black"),
          legend.text  = element_text(family = FAM, size = 9, colour = "black"),
          legend.title = element_text(family = FAM, size = 10, colour = "black"))
  g
}

# Six views: the front view first (the one that matches the 2D panel's orientation),
# then a turn around the vertical axis, then two tilts.
views <- list(
  list(a =   0, e =  0,  t = "front (matches 2D rotate_90)"),
  list(a =  60, e =  0,  t = "azimuth 60°"),
  list(a = 120, e =  0,  t = "azimuth 120°"),
  list(a = 180, e =  0,  t = "azimuth 180°"),
  list(a =  30, e = 30,  t = "azimuth 30°, elevation 30°"),
  list(a =  30, e = -30, t = "azimuth 30°, elevation -30°"))

gs <- lapply(seq_along(views), function(i)
  mk_view(views[[i]]$a, views[[i]]$e, views[[i]]$t, legend = FALSE))

g <- wrap_plots(gs, ncol = 3) +
  plot_annotation(
    title = "3D UMAP, oriented to the rotate_90 convention",
    subtitle = paste("Separate n.components = 3 embedding on integrated.cca dims 1-30;",
                     "not a rotation of the 2D UMAP"),
    theme = theme(plot.title = element_text(family = FAM, face = "bold", size = 13),
                  plot.subtitle = element_text(family = FAM, size = 10)))

# the legend comes from one extra panel so every scatter can stay legend-free
g_leg <- mk_view(0, 0, "front, with key", legend = TRUE)
ggsave(file.path(OUT, "Fig_umap3d_angles.pdf"), g, device = cairo_pdf,
       width = 13.5, height = 9.5)
ggsave(file.path(OUT, "Fig_umap3d_angles.png"), g, width = 13.5, height = 9.5,
       dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_umap3d_front_with_key.pdf"), g_leg, device = cairo_pdf,
       width = 8.5, height = 6.5)
ggsave(file.path(OUT, "Fig_umap3d_front_with_key.png"), g_leg, width = 8.5,
       height = 6.5, dpi = 300, bg = "white")

# each view on its own, for picking one for a figure panel
for (i in seq_along(views)) {
  v <- views[[i]]
  nm <- sprintf("Fig_umap3d_a%03d_e%+03d.png", v$a, v$e)
  ggsave(file.path(OUT, nm), gs[[i]], width = 6.5, height = 6.5, dpi = 300,
         bg = "white")
}
cat("wrote", length(views), "single views + the grid\n")

## ---- optional interactive version -------------------------------------------
if (requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)) {
  p <- plotly::plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~cell_type,
                       colors = PAL_CT, type = "scatter3d", mode = "markers",
                       marker = list(size = 1.6, opacity = 0.75))
  p <- plotly::layout(p, title = "3D UMAP (rotate_90 orientation)",
                      scene = list(xaxis = list(title = "UMAP 1"),
                                   yaxis = list(title = "UMAP 2"),
                                   zaxis = list(title = "UMAP 3")))
  # selfcontained = FALSE deliberately: TRUE needs pandoc, which is not on the
  # compute nodes (job 669544 died there after every figure was already written).
  # The page therefore needs umap3d_interactive_files/ beside it -- copy both, or
  # the plot comes up blank.
  htmlwidgets::saveWidget(p, file.path(OUT, "umap3d_interactive.html"),
                          selfcontained = FALSE)
  cat("wrote interactive HTML (keep umap3d_interactive_files/ next to it)\n")
} else {
  cat("SKIP interactive HTML: plotly and/or htmlwidgets not available.\n")
}

cat("\nwrote ->", OUT, "\n")
cat("DONE\n")
