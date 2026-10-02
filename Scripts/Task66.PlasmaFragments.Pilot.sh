#!/bin/sh
# ================================================================
# Task 66 pilot - rebuild one plasma donor's fragments from raw reads
#   GSM4502069 / SRX8181977 / SRR11615795 (healthy control, GSE149438)
# Date: Oct 1, 2026
#
# Steps (each logged under Data/gse149438_fragment_pilot/, gitignored):
#   0. inputs: FASTQ MD5s match ENA; hg19 from UCSC (hg19.fa.gz, gzip -t)
#   1. trimming (cutadapt 5.2), as recommended for Swift Accel-NGS
#      Methyl-Seq: pass 1 trims adapters (AGATCGGAAGAGC, both reads) and
#      bases below quality 20; pass 2 clips 10 bp from the 5' end of R1,
#      15 bp from the 5' end of R2 and 10 bp from both 3' ends, keeping
#      pairs with both reads at least 30 bp
#   2. Bismark 0.24.2 genome preparation (bowtie2 2.5.5), full hg19
#   3. Bismark paired-end alignment, directional, default scoring
#   4. deduplicate_bismark --paired (position-based, like sambamba's
#      marking in mHapBrowser's pipeline; there are no UMIs)
#   5. fragment extraction, linked and mHap-rule
#      (Task66.PlasmaFragments.Extract.Oct01.2026.py)
#   6. comparison (Task66.PlasmaFragments.Compare.Oct01.2026.R)
# Steps 0 and 1 were run by hand before this script; the commands are
# recorded here.
# ================================================================
set -e
REPO=$HOME/Desktop/bioinformatics-research
P=$REPO/Data/gse149438_fragment_pilot
B=$REPO/Data/shared_fragment_tools/Bismark-0.24.2
cd "$P"

# step 1 (already done):
#   cutadapt -j 2 -q 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 1 --interleaved R1 R2 |
#   cutadapt -j 2 --interleaved -u 10 -U 15 -u -10 -U -10 -m 30 --pair-filter=any -o R1.trim -p R2.trim -

# step 2: wait for genome preparation if it is still running
while pgrep -f bismark_genome_preparation > /dev/null; do sleep 60; done
test -f hg19_bismark/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 || test -f hg19_bismark/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2l
test -f hg19_bismark/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2 || test -f hg19_bismark/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2l
echo "genome ready $(date)"

# step 3
mkdir -p bismark_out tmp
"$B/bismark" --genome hg19_bismark --path_to_bowtie2 /opt/homebrew/bin --samtools_path /opt/homebrew/bin \
  --parallel 1 -p 6 --temp_dir tmp -o bismark_out \
  -1 SRR11615795_1.trim.fq.gz -2 SRR11615795_2.trim.fq.gz > bismark_align.log 2>&1
echo "aligned $(date)"

# step 4
"$B/deduplicate_bismark" --paired --bam --samtools_path /opt/homebrew/bin --output_dir bismark_out \
  bismark_out/SRR11615795_1.trim_bismark_bt2_pe.bam > dedup.log 2>&1
echo "deduplicated $(date)"

# step 5
cd "$REPO"
.venv-shared-fragment/bin/python Scripts/Task66.PlasmaFragments.Extract.Oct01.2026.py \
  "$P/bismark_out/SRR11615795_1.trim_bismark_bt2_pe.deduplicated.bam" "$P/SRR11615795" > "$P/extract.log" 2>&1
echo "extracted $(date)"

# step 6
cd "$REPO/Scripts"
R CMD BATCH --no-save --no-restore Task66.PlasmaFragments.Compare.Oct01.2026.R
echo "compared $(date)"
