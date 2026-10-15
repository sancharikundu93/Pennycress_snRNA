#!/bin/bash
# Step 77 -- top-2 marker dot plot (loads the 5 GB object).
#
#   sbatch run_77_top2_dotplot.sh
#
# Reads one CSV (two GO CSVs (bulk step 31 + snRNA step 03b)) and draws two dotplots. No Seurat, no
# object, no re-enrichment -- it finishes in seconds and is safe to run interactively:
#
#   module reset
#   module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
#   export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2
#   Rscript scripts/77_top2_marker_dotplot.R
#
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=01:00:00
#SBATCH --job-name=dotplot77
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/top2_dotplot_77_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2
unset PYTHONPATH

echo "Host:     $(hostname)"
echo "Start:    $(date)"

Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/77_top2_marker_dotplot.R

echo "End:      $(date)"
