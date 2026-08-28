#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task15.flagbias.noextref
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task15 flag-bias (state + sample 14-17, no external ref), the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save ~/Experiments.Outlier.July31.2026/Scripts/Task15.FlagBias.StateAndSample.NoExtRef.Aug20.2026.R ~/Experiments.Outlier.July31.2026/Scripts/Task15.FlagBias.StateAndSample.NoExtRef.Aug20.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task15 finished, the time is ----------------------------"
date;
