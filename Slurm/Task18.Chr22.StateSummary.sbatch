#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J task18.chr22
#SBATCH --mem=32G

echo "-----------------------------------------------------------------"
echo "Begin running Task18 Chr22 Flag Summary By Methylation State"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task18.Chr22.StateSummary.Aug27.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task18.Chr22.StateSummary.Aug27.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task18 finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
