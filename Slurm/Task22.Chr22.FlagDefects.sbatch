#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J task22.chr22
#SBATCH --mem=32G

echo "-----------------------------------------------------------------"
echo "Begin running Task22 Chr22 Outlier Flag Defect Analysis"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task22.Chr22.FlagDefects.Aug28.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task22.Chr22.FlagDefects.Aug28.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task22 finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
