#!/usr/bin/env Rscript
# Step 21 -- top-N marker dot plot per cell type, split Control vs Salt.
#
# Ports troubleshooting_scripts/scripts_07312026/20_dotplot_top3_res06.R onto the
# nCount > 500 object. That script produced
#   results/dot_plot/res0.6_top3_markers_by_condition.pdf
# from penny_int_res0.6_073126.rds -- the 1000-gate object, 18,475 nuclei, 20
# clusters, no pericycle_like. Construction here is unchanged (split DotPlot + three
# hand-built side legends assembled with patchwork); only the inputs move.
#
# What differs from the old figure, and why:
#   * 23 clusters, not 20 -- top N is taken per CLUSTER then de-duplicated, so the
#     three Endodermis clusters can each contribute their own genes to one row
#   * cluster 5 (pericycle_like) now appears, and Cortex is two rows (Cortex_1,
#     Cortex_2) because annotation curated them as separate cortical layers
#   * markers come from results/02_annotation/markers_RNA_snn_res_0_6.csv, the
#     FindAllMarkers run on this object
#
# TOPN = 2 by default, as in the source script: 3 per cluster gave 60 genes over 36
# rows, denser than the original 42 genes / 28 rows and hard to read. With 23
# clusters here, 2 is the more necessary choice, not less.
#
#   sbatch run_21_dotplot_top2.sh
#   TOPN=3 sbatch --export=ALL run_21_dotplot_top2.sh
#
# Outputs -> ncount500_08132026/results/21_dotplot_top_markers/
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
options(bitmapType = "cairo")

FAM   <- "Times New Roman"
INK   <- "black"
C_COL <- "#F0BC35"                  # Control
S_COL <- "#C2185B"                  # Salt

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
RDS  <- file.path(ROOT, "results/penny_mn106_ncount500_integrated.rds")
MK   <- file.path(ROOT, "results/02_annotation/markers_RNA_snn_res_0_6.csv")
LAB  <- file.path(ROOT, "results/02_annotation/labels_res0.6.csv")
OUT  <- file.path(ROOT, "results/21_dotplot_top_markers")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/palette_celltypes.R"))
stopifnot(file.exists(RDS), file.exists(MK), file.exists(LAB))

TOPN <- as.integer(Sys.getenv("TOPN", "2"))

L   <- read.csv(LAB, colClasses = "character")
lab <- setNames(L$cell_type, L$cluster)
# CT_LEVELS is the shared anatomical order (outer -> stele), so this figure's rows
# read in the same direction as the proportion bar and the ridge plots.
ct_levels <- CT_LEVELS[CT_LEVELS %in% unique(lab)]
stopifnot(setequal(ct_levels, unique(lab)))

so <- readRDS(RDS)
cl <- as.character(so@meta.data[["RNA_snn_res.0.6"]]); names(cl) <- colnames(so)
if (any(is.na(lab[cl]))) stop("unlabelled cluster(s): ",
                              paste(sort(unique(cl[is.na(lab[cl])])), collapse = ", "))
ct <- factor(lab[cl], levels = rev(ct_levels)); names(ct) <- colnames(so)
so$cluster_label <- ct
so$condition <- factor(as.character(so$condition), levels = c("Control", "Salt"))
Idents(so) <- "cluster_label"
cat("nuclei:", ncol(so), " cell types:", nlevels(droplevels(so$cluster_label)), "\n")

markers <- read.csv(MK)
topn <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = TOPN) %>% ungroup()

# Order the x axis to follow the y axis so the marker blocks form one clean diagonal.
# y levels are rev(ct_levels) and ggplot puts the FIRST level at the bottom, so genes
# grouped bottom-to-top give a diagonal running bottom-left to top-right.
ct_order <- rev(ct_levels)
topn <- topn %>%
  mutate(celltype = lab[as.character(cluster)],
         ct_rank  = match(celltype, ct_order)) %>%
  arrange(ct_rank, as.integer(as.character(cluster)), desc(avg_log2FC))
top_genes <- unique(topn$gene)
top_genes <- top_genes[top_genes %in% rownames(so)]
top_clean <- gsub("\\.v4\\.1$", "", top_genes)
names(top_clean) <- top_genes
cat("marker genes on the x axis:", length(top_genes),
    "(from", nrow(topn), "cluster-gene pairs;",
    nrow(topn) - length(top_genes), "shared between clusters)\n")
write.csv(topn[, c("cluster", "celltype", "gene", "avg_log2FC", "pct.1", "pct.2",
                   "p_val_adj")],
          file.path(OUT, sprintf("top%d_markers_per_cluster.csv", TOPN)),
          row.names = FALSE)

combined_dot <- DotPlot(so,
                        features = top_genes,
                        group.by = "cluster_label",
                        split.by = "condition",
                        cols = c("Control" = C_COL, "Salt" = S_COL),
                        col.min = 0,
                        dot.scale = 10,
                        scale = TRUE) +
  scale_x_discrete(labels = top_clean) +
  labs(title = sprintf("Top %d Marker Genes per Cell Type", TOPN),
       x = "Gene", y = "Cell Type") +
  theme(text = element_text(family = FAM, colour = INK),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 19,
                                   family = FAM, colour = INK),
        axis.text.y = element_text(size = 20, family = FAM, colour = INK),
        axis.title = element_text(size = 21, family = FAM, colour = INK),
        plot.title = element_text(size = 25, face = "bold", family = FAM, colour = INK),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position = "none",
        panel.spacing = unit(0, "lines"))

# ---- y axis: one centred name per Control/Salt pair --------------------------
# DotPlot(split.by=) emits "<cell type>_Control" and "<cell type>_Salt", so every
# name is printed twice with a suffix. Blank the axis text and draw one name at each
# pair midpoint outside the panel instead -- a discrete scale can only put a label
# ON a row, never between two.
ylev  <- levels(combined_dot$data$id)
ct_of <- sub("_(Control|Salt)$", "", ylev)
mid   <- tapply(seq_along(ylev), factor(ct_of, levels = unique(ct_of)), mean)
ylab_df <- data.frame(y = as.numeric(mid), lab = names(mid))

combined_dot <- combined_dot +
  geom_text(data = ylab_df, aes(x = 0.35, y = y, label = lab),
            inherit.aes = FALSE, hjust = 1, size = 7.6,
            family = FAM, colour = INK) +
  coord_cartesian(clip = "off") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # the names are drawn with clip = "off", so ggplot reserves no room for them
        # and "Cell Type" collides; reserve the gap as the title's right margin
        axis.title.y = element_text(margin = margin(r = 215)),
        plot.margin = margin(6, 6, 6, 8))

# ---- side legends -------------------------------------------------------------
# Seurat's split DotPlot maps each condition from GREY to that condition's colour --
# there is no white midpoint, so these are plain two-stop ramps that match what the
# dots actually do.
grad <- function(high, title) {
  # geom_tile with an explicit width, NOT geom_raster: raster stretches to the panel
  # width and draws the ramp as a block with the border floating inside it.
  ggplot(data.frame(y = seq(0, 1, length.out = 256)), aes(x = 1, y = y, fill = y)) +
    geom_tile(width = 0.30, height = 1 / 255) +
    scale_fill_gradient(low = "grey85", high = high) +
    annotate("rect", xmin = 0.85, xmax = 1.15, ymin = 0, ymax = 1,
             fill = NA, colour = INK, linewidth = 0.4) +
    annotate("text", x = 1.26, y = 1, label = "Max", hjust = 0, vjust = 0.5,
             size = 4.2, family = FAM, colour = INK) +
    annotate("text", x = 1.26, y = 0, label = "Min", hjust = 0, vjust = 0.5,
             size = 4.2, family = FAM, colour = INK) +
    annotate("text", x = 0.85, y = 1.14, label = title, hjust = 0, vjust = 0,
             size = 5.0, fontface = "bold", family = FAM, colour = INK) +
    coord_cartesian(xlim = c(0.80, 2.45), ylim = c(-0.06, 1.34), clip = "off") +
    theme_void(base_family = FAM) +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
}
control_grad <- grad(C_COL, "Control")
salt_grad    <- grad(S_COL, "Salt")

# Dot size key on the identical scale DotPlot uses (scale_radius over pct expressed),
# so the keys are the same size as the dots they explain.
pct_vals <- c(0, 10, 20, 30, 40, 50)
pct_legend <- ggplot(data.frame(size = pct_vals, y = seq_along(pct_vals), x = 1),
                     aes(x, y, size = size)) +
  geom_point(colour = INK) +
  scale_radius(limits = c(0, 100), range = c(0, 10)) +
  annotate("text", x = 1.55, y = seq_along(pct_vals), label = pct_vals,
           hjust = 0, size = 4.6, family = FAM, colour = INK) +
  annotate("text", x = 0.55, y = length(pct_vals) + 1.1, label = "Pct expressed",
           hjust = 0, vjust = 0, size = 5.4, fontface = "bold",
           family = FAM, colour = INK) +
  coord_cartesian(xlim = c(0.55, 2.6), ylim = c(0.2, length(pct_vals) + 1.9),
                  clip = "off") +
  theme_void(base_family = FAM) +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))

all_legends <- plot_spacer() / control_grad / salt_grad / pct_legend / plot_spacer() +
  plot_layout(heights = c(1.0, 1.7, 1.7, 2.1, 1.0))

final_plot <- combined_dot + all_legends + plot_layout(widths = c(17, 0.95))

# ~0.5 in per gene column; each cell type occupies two rows (Control + Salt) and the
# pair only has to be tall enough for the dots.
W <- max(20, 0.50 * length(top_genes) + 6)
H <- max(11, 0.29 * (2 * nlevels(droplevels(so$cluster_label))) + 2)
stem <- file.path(OUT, sprintf("Fig_dotplot_top%d_by_condition", TOPN))
ggsave(paste0(stem, ".pdf"), final_plot, width = W, height = H,
       device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(stem, ".png"), final_plot, width = W, height = H, dpi = 200,
       limitsize = FALSE)
cat("wrote to", OUT, " width =", round(W, 1), " height =", round(H, 1), "\n")
cat("DONE\n")
