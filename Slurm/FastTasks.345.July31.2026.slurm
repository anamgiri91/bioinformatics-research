#!/bin/bash
#SBATCH --job-name=fast_tasks_345
#SBATCH --partition=shared
#SBATCH -n 1
#SBATCH --output=/home/wln26/fast_tasks_345_%j.log
#SBATCH --error=/home/wln26/fast_tasks_345_%j.log
#SBATCH --time=00:20:00
#SBATCH --mem=8G

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files
Rscript ~/fast_tasks_3_4_5.R
