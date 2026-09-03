#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -J task33.genomewide
#SBATCH --mem=64G

echo "-----------------------------------------------------------------"
echo "Begin running Task33 genome-wide headroom, geometry, stability, burden"
echo "Start time:"
date
echo "-----------------------------------------------------------------"

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files

R CMD BATCH --no-save --no-restore \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task33.GenomeWideHeadroom.Sep03.2026.R \
/home/wln26/Experiments.Outlier.July31.2026/Scripts/Task33.GenomeWideHeadroom.Sep03.2026.Rout

echo "-----------------------------------------------------------------"
echo "task33.genomewide finished"
echo "End time:"
date
echo "-----------------------------------------------------------------"
