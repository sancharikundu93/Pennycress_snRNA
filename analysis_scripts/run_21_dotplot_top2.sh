#!/bin/bash
# Step 21 -- top-2 marker dot plot per cell type, split Control vs Salt.
#   sbatch run_21_dotplot_top2.sh
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --job-name=dot21
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/dotplot_top_%j.log

set -eo pipefail
mkdir -p /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs

module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH

echo "Host: $(hostname)"
Rscript /projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/scripts/21_dotplot_top2_by_condition.R

echo "21_dotplot_top2 done"
