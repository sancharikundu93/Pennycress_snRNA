#!/bin/bash
# Step 1 -- export the slim app data from the 5 GB integrated object.
#   sbatch run_step1_export_app_data.sh
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --job-name=app_data
#SBATCH --output=/projects/songli_lab/A_Sanchari/A_pennycress_snRNA/ncount500_08132026/results/logs/app_data_%j.log

set -eo pipefail
module reset
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load GLPK/5.0-GCCcore-13.3.0
unset PYTHONPATH
export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2

echo "Host: $(hostname)"
Rscript /home/sancharikundu/github/pennycress_snrna/data_prep/step1_export_app_data.R
