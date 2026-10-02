#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name=task68.sperm.qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-2%2
#SBATCH --output=task68-sperm-%A_%a.log
set -euo pipefail
# Activate a verified environment before sbatch (Python>=3.11, numpy, pysam,
# Trim Galore, cutadapt, Bismark, Bowtie2 and samtools). No plasma processing.
: "${SRN_REPO:?set repository root}"
: "${SRN_SCRATCH:?set private scratch directory}"
: "${SRN_REFERENCE_DIR:?set full hg38 Bismark reference directory}"
: "${SRN_CPG_MAP:?set CpG map built from the same reference FASTA}"
cd "$SRN_REPO"
python3 Scripts/Task68.SpermPilot.py \
  --manifest Results/Task68_SpermPilotSamples.tsv --row "$SLURM_ARRAY_TASK_ID" \
  --work "$SRN_SCRATCH" --reference-dir "$SRN_REFERENCE_DIR" --map-dir "$SRN_CPG_MAP" \
  --threads "$SLURM_CPUS_PER_TASK"
