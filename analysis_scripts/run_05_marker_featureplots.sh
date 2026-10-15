#!/bin/bash
# Step 05 -- canonical marker FeaturePlots (Control | Salt) + marker reference table.
#   sbatch run_05_marker_featureplots.sh                              # literature markers
#   MARKER_SET=all sbatch --export=ALL run_05_marker_featureplots.sh  # full 45-marker panel
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --job-name=marker_fplots
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/marker_featureplots_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
echo "MARKER_SET=${MARKER_SET:-literature}"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/05_canonical_marker_featureplots.R

echo "05_marker_featureplots done"
