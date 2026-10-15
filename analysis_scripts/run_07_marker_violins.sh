#!/bin/bash
# Step 07 -- canonical marker violins + dot plot, grouped by cell type (and cluster).
#   sbatch run_07_marker_violins.sh
#   MARKER_SET=all sbatch --export=ALL run_07_marker_violins.sh
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --job-name=marker_violins
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/marker_violins_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
echo "MARKER_SET=${MARKER_SET:-literature}"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/07_canonical_marker_violins.R

echo "07_marker_violins done"
