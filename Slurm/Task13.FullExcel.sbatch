#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 120:00:00
#SBATCH -J task13.fullexcel
#SBATCH --mem=32G

echo "-----------------------------------------------------------------"
echo "Begin running Task13 full Excel matrix generation"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save ~/Experiments.Outlier.July31.2026/Scripts/Task13.CreateFullExcel.R ~/Experiments.Outlier.July31.2026/Scripts/Task13.CreateFullExcel.Rout

echo "-----------------------------------------------------------------"
echo "Task13 full Excel matrix finished"
date;
