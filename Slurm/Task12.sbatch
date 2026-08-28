#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task12.external.ref
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task12 external reference comparison, the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH ~/Experiments.Outlier.July31.2026/Scripts/Task12.ExternalReference.Aug6.2026.R ~/Experiments.Outlier.July31.2026/Scripts/Task12.ExternalReference.Aug6.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task12 finished, the time is ----------------------------"
date;
