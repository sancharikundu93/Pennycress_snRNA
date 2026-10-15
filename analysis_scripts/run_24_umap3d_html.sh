#!/bin/bash
# Step 24 -- self-contained rotating 3D UMAP HTML.
#   sbatch run_24_umap3d_html.sh
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --job-name=html24
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/umap3d_html_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/24_umap3d_interactive_html.R

echo "24_umap3d_html done"
