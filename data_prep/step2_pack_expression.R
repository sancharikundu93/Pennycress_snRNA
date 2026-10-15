#!/usr/bin/env Rscript
# Step 2 -- pack the expression matrix so the app reads one gene at a time.
#
#   Rscript step2_pack_expression.R
#
# Held in memory, the cells x genes matrix (48M nonzeros) costs ~600 MB, most of the
# 1 GB a free shinyapps.io instance gets. So it goes to disk instead, one block per gene:
#   int32 cell indices (1-based) then int16 values (log-normalised x 1000)
# and the app seeks straight to the gene it needs. Values keep 3 decimals, which is
# what step 1 rounded to already; the int16 ceiling (32.767) is far above any
# log-normalised value.
#
# In:  results/84_shiny_app_data/app_data.rds   (from step 1)
# Out: app/data/meta.rds, app/data/expr.bin
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
IN  <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/84_shiny_app_data/app_data.rds"
OUT <- "/home/sancharikundu/github/pennycress_snrna/app/data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

D <- readRDS(IN)
X <- D$expr                                   # cells x genes dgCMatrix
stopifnot(max(X@x) < 32.767)
con <- file(file.path(OUT, "expr.bin"), "wb")
n   <- diff(X@p)
off <- numeric(ncol(X))
pos <- 0
for (j in seq_len(ncol(X))) {
  k <- if (n[j]) (X@p[j] + 1):X@p[j + 1] else integer(0)
  off[j] <- pos
  writeBin(as.integer(X@i[k] + 1L), con, size = 4)
  writeBin(as.integer(round(X@x[k] * 1000)), con, size = 2)
  pos <- pos + n[j] * 6
}
close(con)
D$genes[, `:=`(offset = off, nnz = n)]
D$expr <- NULL
D$n_cells <- nrow(X)
saveRDS(D, file.path(OUT, "meta.rds"))
cat(sprintf("genes %d, expr.bin %.0f MB, meta.rds %.1f MB\n", ncol(X),
            file.size(file.path(OUT, "expr.bin")) / 1e6, file.size(file.path(OUT, "meta.rds")) / 1e6))
