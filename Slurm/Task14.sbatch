#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task14.flagbias.extRef
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task14 flag-bias (state + sample14-17, external ref), the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save ~/Experiments.Outlier.July31.2026/Scripts/Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R ~/Experiments.Outlier.July31.2026/Scripts/Task14.FlagBias.StateAndSample.ExtRef.Aug2026.Rout

echo "-----------------------------------------------------------------"
echo "Task14 finished, the time is ----------------------------"
date;
