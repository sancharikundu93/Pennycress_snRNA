#!/bin/bash
# Step 23 -- 3D UMAP oriented to the rotate_90 convention.
#   sbatch run_23_umap3d.sh
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=4
#SBATCH --mem=250G
#SBATCH --time=03:00:00
#SBATCH --job-name=umap3d23
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/umap3d_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/23_umap3d_rotate90.R

echo "23_umap3d done"
