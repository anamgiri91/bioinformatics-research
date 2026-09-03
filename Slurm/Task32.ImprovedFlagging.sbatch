#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J task32.improve
#SBATCH --mem=48G

echo "-----------------------------------------------------------------"
echo "Begin running Task32 measured noise and improved flagging rules"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save --no-restore \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task32.ImprovedFlagging.Sep03.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task32.ImprovedFlagging.Sep03.2026.Rout

echo "-----------------------------------------------------------------"
echo "task32.improve finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
