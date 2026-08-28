#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task17.chr22.validation
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task17 chr22 biological/local CpG validation"
echo "The time is:"
date

echo "-----------------------------------------------------------------"
echo "Working directory:"
cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files
pwd

echo "-----------------------------------------------------------------"
echo "Running Task17 R script..."
echo "-----------------------------------------------------------------"

R CMD BATCH --no-save \
~/Experiments.Outlier.July31.2026/Scripts/Task17.Chr22.BiologicalValidation.Aug21.2026.R \
~/Experiments.Outlier.July31.2026/Scripts/Task17.Chr22.BiologicalValidation.Aug21.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task17 chr22 validation finished"
echo "The time is:"
date

echo "-----------------------------------------------------------------"
echo "Job finished."
echo "-----------------------------------------------------------------"
