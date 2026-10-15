#!/usr/bin/env Rscript
# Step 1 -- export the slim data file the Shiny explorer reads.
#
#   sbatch run_step1_export_app_data.sh
#
# The integrated Seurat object is 5 GB on disk and needs ~20 GB to load; the free
# shinyapps.io instance has 1 GB of RAM. So the app never sees the object. This
# writes one .rds holding only what the explorer draws:
#   cells  -- barcode, sample, condition, cluster, cell type, 2D + 3D UMAP coords
#   genes  -- feature id, locus, Arabidopsis best-hit AGI, TAIR symbol(s)
#   expr   -- log-normalised RNA counts, CELLS x GENES sparse, so pulling one gene
#             is a single column slice (fast) rather than a row slice (slow)
#
# Cell-type labels are the step-74 from-scratch annotation (clusters 15 and 22
# relabelled), not the older labels baked into umap3d_embedding.csv.
#
# Outputs -> ncount500_08132026/results/84_shiny_app_data/app_data.rds
.libPaths(c("/home/sancharikundu/R/owl-genoa/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(org.At.tair.db); library(AnnotationDbi)
})

LIN  <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026"
RDS  <- file.path(LIN, "results/penny_mn106_ncount500_integrated.rds")
LAB  <- file.path(LIN, "results/74_annotation_from_scratch/labels_res0.6_from_scratch.csv")
U3D  <- file.path(LIN, "results/23_umap3d/umap3d_embedding.csv")
MAP  <- file.path(LIN, "results/81_agi_locus_map/agi_to_locus_map.csv.gz")
OUT  <- file.path(LIN, "results/84_shiny_app_data")
source(file.path(LIN, "scripts/palette_celltypes.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c(RDS, LAB, U3D, MAP)) stopifnot(file.exists(f))

cat("=== step 1: export app data ===\n")
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
cat("object:", ncol(obj), "nuclei x", nrow(obj), "genes\n")
cat("meta columns:", paste(colnames(obj@meta.data), collapse = ", "), "\n")
cat("reductions:", paste(Reductions(obj), collapse = ", "), "\n")

## ---- cells --------------------------------------------------------------------
cl  <- as.character(obj$RNA_snn_res.0.6)
lab <- fread(LAB)
hit <- match(cl, as.character(lab$cluster))
if (anyNA(hit)) stop("unlabelled cluster(s): ", paste(unique(cl[is.na(hit)]), collapse = ","))
ct  <- lab$cell_type[hit]

u2 <- Embeddings(obj, "umap.integrated")
u3 <- fread(U3D)
# umap3d_embedding.csv has no barcode column; step 23 wrote it in object order.
# Verify that rather than trust it: every row's cluster must agree.
stopifnot(nrow(u3) == ncol(obj), all(as.character(u3$cluster) == cl))

samp <- as.character(obj$sample)
cells <- data.table(cell = colnames(obj), sample = samp,
                    condition = as.character(obj$condition), cluster = cl,
                    cell_type = ct,
                    u1 = u2[, 1], u2 = u2[, 2],
                    x3 = u3$x, y3 = u3$y, z3 = u3$z)
CT_ORD <- CT_LEVELS[CT_LEVELS %in% unique(ct)]
CT_ORD <- c(CT_ORD, setdiff(sort(unique(ct)), CT_ORD))
cat("cell types:", length(CT_ORD), "\n"); print(cells[, .N, by = .(condition)])

## ---- expression ---------------------------------------------------------------
if (length(Layers(obj[["RNA"]], search = "data")) > 1) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
X <- LayerData(obj, assay = "RNA", layer = "data")
X <- as(X, "CsparseMatrix")
keep <- Matrix::rowSums(X > 0) >= 3          # drop genes seen in fewer than 3 nuclei
X <- X[keep, , drop = FALSE]
X@x <- round(X@x, 3)
expr <- t(X)                                 # cells x genes, dgCMatrix
cat("genes kept (>=3 nuclei):", ncol(expr), " nonzeros:", length(expr@x), "\n")
rm(obj, X); invisible(gc())

## ---- genes --------------------------------------------------------------------
feat  <- colnames(expr)
locus <- sub("\\.v4\\.1$", "", feat)
mp    <- fread(MAP)[route == "best_hit"]
agi   <- mp[, .(agi = paste(sort(unique(agi)), collapse = ";")), by = locus]
genes <- data.table(feature = feat, locus = locus)
genes <- merge(genes, agi, by = "locus", all.x = TRUE, sort = FALSE)
genes <- genes[match(feat, genes$feature)]

all_agi <- unique(unlist(strsplit(na.omit(genes$agi), ";")))
sym <- as.data.table(AnnotationDbi::select(org.At.tair.db, keys = all_agi,
                                           keytype = "TAIR", columns = "SYMBOL"))
sym <- sym[!is.na(SYMBOL)]
sym_by_agi <- sym[, .(symbols = paste(unique(SYMBOL), collapse = ";")), by = TAIR]
first_sym  <- function(a) {
  if (is.na(a)) return(NA_character_)
  s <- sym_by_agi$symbols[match(strsplit(a, ";")[[1]], sym_by_agi$TAIR)]
  s <- s[!is.na(s)]
  if (!length(s)) NA_character_ else paste(unique(unlist(strsplit(s, ";"))), collapse = ";")
}
genes[, symbols := vapply(agi, first_sym, character(1))]
genes[, symbol := sub(";.*", "", symbols)]
cat("genes with AGI:", sum(!is.na(genes$agi)), " with symbol:", sum(!is.na(genes$symbol)), "\n")

app <- list(cells = cells, genes = genes, expr = expr, ct_levels = CT_ORD,
            ct_colours = ct_palette(CT_ORD),
            made = format(Sys.time(), "%Y-%m-%d"))
saveRDS(app, file.path(OUT, "app_data.rds"), compress = "xz")
cat("wrote", file.path(OUT, "app_data.rds"),
    sprintf("(%.0f MB)\n", file.size(file.path(OUT, "app_data.rds")) / 1e6))
