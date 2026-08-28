#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J task17.chr22
#SBATCH --mem=16G

echo "-----------------------------------------------------------------"
echo "Begin running Task17 Chr22 Biological Validation"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task17.Chr22.BiologicalValidation.Aug21.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task17.Chr22.BiologicalValidation.Aug21.2026.Rout

echo "-----------------------------------------------------------------"
echo "Task17 finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
