#!/bin/bash
set -euo pipefail
: "${SRN_REPO:?set repository root}"
: "${SRN_SCRATCH:?set private scratch directory}"
: "${SRN_REFERENCE_DIR:?set full hg38 Bismark reference directory}"
: "${SRN_CPG_MAP:?set CpG map directory}"
for tool in python3 curl trim_galore bismark deduplicate_bismark bowtie2 samtools sbatch; do
  command -v "$tool"
done
python3 -c 'import sys,numpy,pysam; assert sys.version_info >= (3,11); print(sys.version); print("numpy",numpy.__version__,"pysam",pysam.__version__)'
samtools --version
bowtie2 --version
bismark --version
trim_galore --version
mkdir -p "$SRN_SCRATCH"
df -h "$SRN_SCRATCH"
test -f "$SRN_CPG_MAP/manifest.json"
cd "$SRN_REPO"
python3 Scripts/Task68.SpermPilot.py --manifest Results/Task68_SpermPilotSamples.tsv --row 1 \
  --work "$SRN_SCRATCH" --reference-dir "$SRN_REFERENCE_DIR" --map-dir "$SRN_CPG_MAP" --dry-run
