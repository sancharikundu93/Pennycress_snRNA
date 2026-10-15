#!/bin/bash
# Step 78 -- thresholded canonical-marker feature plot grid.
#
#   sbatch run_78_featureplot_grid.sh
#
# Reads one CSV (two GO CSVs (bulk step 31 + snRNA step 03b)) and draws two dotplots. No Seurat, no
# object, no re-enrichment -- it finishes in seconds and is safe to run interactively:
#
#   module reset
#   module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
#   export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2
#   Rscript scripts/78_marker_featureplot_grid.R
#
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=01:00:00
#SBATCH --job-name=feat6
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/featureplot_6_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2
unset PYTHONPATH

echo "Host:     $(hostname)"
echo "Start:    $(date)"

export GENES="CASP1,ESB1,PRP3,AGP3,APL,XCP1"
export NCOL=4
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/78_marker_featureplot_grid.R

echo "End:      $(date)"
