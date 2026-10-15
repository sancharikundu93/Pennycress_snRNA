#!/bin/bash
# Step 02 -- markers + annotation evidence + cross-tab against the 1000-gate labels.
#   CLUSTER_COL=RNA_snn_res.0.6 sbatch --export=ALL run_02_markers_annotate.sh
#
# Memory: loads the new 3.5 GB object, runs FindAllMarkers, then loads the OLD 3.7 GB
# object for the cross-tab. They are read sequentially with an rm()+gc() between, so
# peak is one object plus the marker tables -- 250G is comfortable.
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=08:00:00
#SBATCH --job-name=annot_ncount500
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/markers_annotate_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
echo "CLUSTER_COL=${CLUSTER_COL:-RNA_snn_res.0.6}"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/02_markers_annotate_ncount500.R

echo "02_markers_annotate done"
