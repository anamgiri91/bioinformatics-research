#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task10.flag.state.table
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task10 flag-by-state table script, the time is ----------------------------"
date;

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH ~/Experiments.Outlier.July31.2026/Scripts/Task10.FlagStateTable.Aug6.2026.R ~/Experiments.Outlier.July31.2026/Scripts/Task10.FlagStateTable.Aug6.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task10 finished, the time is ----------------------------"
date;
