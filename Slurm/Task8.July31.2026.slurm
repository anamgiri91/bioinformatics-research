#!/bin/bash
#SBATCH --job-name=task8_filter
#SBATCH --partition=shared
#SBATCH -n 1
#SBATCH --output=/home/wln26/task8_%j.log
#SBATCH --error=/home/wln26/task8_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G

cd /home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files
Rscript ~/task8_filter.R
