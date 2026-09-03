#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J task29.contrary
#SBATCH --mem=32G

echo "-----------------------------------------------------------------"
echo "Begin running Task29 contrary-flag decomposition over sites and samples"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save --no-restore \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task29.ContraryFlagDecomposition.Sep03.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task29.ContraryFlagDecomposition.Sep03.2026.Rout

echo "-----------------------------------------------------------------"
echo "task29.contrary finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
