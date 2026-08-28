#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -J task9.boxplot.hist
#SBATCH --mem=8G

echo "-----------------------------------------------------------------"
echo "Begin running Task9 boxplot/histogram script, the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH ~/Experiments.Outlier.July31.2026/Scripts/Task9.BoxplotHistogram.Aug6.2026.R ~/Experiments.Outlier.July31.2026/Scripts/Task9.BoxplotHistogram.Aug6.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task9 finished, the time is ----------------------------"
date;
