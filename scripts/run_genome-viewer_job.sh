#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --mem=16G
#SBATCH -J genome-viewer
#SBATCH -p all
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%j.out

# load R module on cluster
module load R/4.2.1

# run the R script
Rscript scripts/genome-viewer.R EGFR --files="/cluster/projects/wouterslab/RW555/deeptools/heatmap/RW555-1_treat_pileup.bw,/cluster/projects/wouterslab/RW555/deeptools/heatmap/RW555-2_treat_pileup.bw"
