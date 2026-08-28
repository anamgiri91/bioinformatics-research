#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task13.sampleflags
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task13 sample-level flag summary, the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save ~/Experiments.Outlier.July31.2026/Scripts/Task13.CreateSampleSummary.R ~/Experiments.Outlier.July31.2026/Scripts/Task13.CreateSampleSummary.Rout

echo "-----------------------------------------------------------------"
echo "Task13 sample-level flag summary finished, the time is ----------------------------"
date;
