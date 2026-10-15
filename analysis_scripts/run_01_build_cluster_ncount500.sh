#!/bin/bash
# Step 07.1 -- rebuild the integrated object at the nCount 500 gate.
#   sbatch run_01_build_cluster_ncount500.sh
#
# Memory: 25,206 nuclei x 27,768 genes. ScaleData runs on ALL genes and densifies to
# ~5.6 GB; with working copies and CCA the peak sits well under 200G. step_06 used
# 150G for the 18,475-nucleus object; 200G here for the extra 6,731 nuclei. Do not
# inflate this further -- an oversized reservation would not schedule when the
# cluster was busy (see the note in step_06/run_01_build_cluster.sh).
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=build_ncount500
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/build_cluster_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
# igraph, via Seurat, needs libglpk.so.40 -- without this the job dies at library(Seurat)
module load GLPK/5.0-GCCcore-13.3.0
# the R module puts its own python on PATH; unset PYTHONPATH so reticulate-backed
# UMAP does not pick up an incompatible environment
unset PYTHONPATH

echo "Host: $(hostname)"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/01_build_cluster_ncount500.R

echo "01_build_cluster_ncount500 done"
