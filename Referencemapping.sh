#!/bin/bash
#$ -N Referencemapping
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=200:00:00
#$ -l s_vmem=8G,mem_req=8G
#$ -pe def_slot 4
#$ -o logs/Referencemapping.out
#$ -e logs/Referencemapping.err

set -e
set -x
set -o pipefail

# ============================================================
# CONFIG
# ============================================================

MAP_ROOT="${MAP_ROOT:-$PWD/results/reference_mapping}"
SAMPLES_ROOT="${SAMPLES_ROOT:-$PWD/data/fastq_by_sample}"
LOG_ROOT="$MAP_ROOT/logs"

REF_SOURCE="${REF_SOURCE:-$PWD/references/reference_genome.fna}"
[ -f "$REF_SOURCE" ] || { echo "ERROR: reference not found: $REF_SOURCE"; exit 1; }

REF_DIR="$MAP_ROOT/reference"
REF_FA="$REF_DIR/$(basename "$REF_SOURCE")"

SIF_DIR="${SIF_DIR:-$HOME/sif_images}"
BWA_SIF="$SIF_DIR/bwa.sif"
SAMTOOLS_SIF="$SIF_DIR/samtools.sif"
MOSDEPTH_SIF="$SIF_DIR/mosdepth.sif"
PICARD_JAR="${PICARD_JAR:-$HOME/tools/picard.jar}"

JAVA_OPTS="-Xms1g -Xmx3g -XX:CompressedClassSpaceSize=128m -XX:MaxMetaspaceSize=512m -Djava.util.concurrent.ForkJoinPool.common.parallelism=4"

SCRATCH_BASE="${TMPDIR:-$MAP_ROOT/tmp_scratch}"

# ------------------------------------------------------------
# Normalization settings
# REFERENCE_MASTER: master_refmap_metrics.csv from the paired
# dataset to normalise against (set to the shallower dataset).
# ------------------------------------------------------------

MIN_MAPPED_PCT=95
MIN_MAPPED_READS=200000
NORMALIZATION_SEED=42
REFERENCE_MASTER="${REFERENCE_MASTER:-$PWD/reference_mapping_comparator/master_refmap_metrics.csv}"
[ -f "$REFERENCE_MASTER" ] || { echo "ERROR: REFERENCE_MASTER not found: $REFERENCE_MASTER"; exit 1; }

# ============================================================
# ENV / MODULES
# ============================================================

module use /usr/local/package/modulefiles
module load apptainer

THREADS="${NSLOTS:-8}"

export _JAVA_OPTIONS="$JAVA_OPTS"
export JAVA_FORKJOINPOOL_COMMON_PARALLELISM=4
ulimit -u 4096 || true

mkdir -p "$MAP_ROOT" "$REF_DIR" "$SCRATCH_BASE" "$LOG_ROOT"

# ============================================================
# STEP -1: Stage reference locally + build indexes (BWA + FAI)
# ============================================================

[ -f "$BWA_SIF" ]      || { echo "ERROR: bwa SIF not found: $BWA_SIF"; exit 1; }
[ -f "$SAMTOOLS_SIF" ] || { echo "ERROR: samtools SIF not found: $SAMTOOLS_SIF"; exit 1; }
[ -f "$MOSDEPTH_SIF" ] || { echo "ERROR: mosdepth SIF not found: $MOSDEPTH_SIF"; exit 1; }
[ -f "$PICARD_JAR" ]   || { echo "ERROR: picard jar not found: $PICARD_JAR"; exit 1; }

if [ ! -f "$REF_FA" ]; then
  cp -f "$REF_SOURCE" "$REF_FA"
fi

for ext in amb ann bwt pac sa; do
  if [ -f "${REF_SOURCE}.${ext}" ] && [ ! -f "${REF_FA}.${ext}" ]; then
    cp -f "${REF_SOURCE}.${ext}" "${REF_FA}.${ext}"
  fi
done

for ext in 0123 0123.bwt 0123.pac 0123.ann 0123.amb 0123.sa; do
  if [ -f "${REF_SOURCE}.${ext}" ] && [ ! -f "${REF_FA}.${ext}" ]; then
    cp -f "${REF_SOURCE}.${ext}" "${REF_FA}.${ext}"
  fi
done

if [ ! -f "${REF_FA}.fai" ]; then
  echo "INFO: Building samtools faidx for $REF_FA"
  singularity exec "$SAMTOOLS_SIF" samtools faidx "$REF_FA"
fi

if [ ! -f "${REF_FA}.bwt" ] && [ ! -f "${REF_FA}.0123" ]; then
  echo "INFO: Building bwa index for $REF_FA"
  singularity exec "$BWA_SIF" bwa index "$REF_FA"
fi

echo "INFO: Reference staging/indexing done. REF_FA=$REF_FA"

[ -f "$REF_FA" ] || { echo "ERROR: project reference copy missing: $REF_FA"; exit 1; }
if [ ! -f "${REF_FA}.bwt" ] && [ ! -f "${REF_FA}.0123" ]; then
  echo "ERROR: BWA index files not found next to reference"; exit 1
fi

REF_FAI="${REF_FA}.fai"
REF_DICT="${REF_FA%.fna}.dict"

if [ ! -f "$REF_FAI" ]; then
  singularity exec "$SAMTOOLS_SIF" samtools faidx "$REF_FA"
fi

if [ ! -f "$REF_DICT" ]; then
  java -jar "$PICARD_JAR" CreateSequenceDictionary R="$REF_FA" O="$REF_DICT"
fi

BIN_10KB="$MAP_ROOT/bins_10kb.bed"
BIN_50KB="$MAP_ROOT/bins_50kb.bed"

if [ ! -s "$BIN_10KB" ]; then
  awk 'BEGIN{OFS="\t"} {for(i=0;i<$2;i+=10000){e=i+10000; if(e>$2)e=$2; print $1,i,e}}' \
    "$REF_FAI" > "$BIN_10KB"
fi

if [ ! -s "$BIN_50KB" ]; then
  awk 'BEGIN{OFS="\t"} {for(i=0;i<$2;i+=50000){e=i+50000; if(e>$2)e=$2; print $1,i,e}}' \
    "$REF_FAI" > "$BIN_50KB"
fi

# ============================================================
# STEP 0: Initialise master CSV headers (once, if not yet created)
# ============================================================

MASTER_REFMAP_CSV="$MAP_ROOT/master_refmap_metrics.csv"
if [ ! -f "$MASTER_REFMAP_CSV" ]; then
  echo "Sample,MAPPED_PCT,PROPERLY_PAIRED_PCT,FRAC_COVERED_GT0,FRAC_COVERED_GE5,FRAC_COVERED_GE10,FRAC_ZERO_DEPTH,MEAN_DEPTH_ALL_POS,MEAN_DEPTH_COVERED_POS,DUP_READS,TOTAL_READS,MAPPED_READS" > "$MASTER_REFMAP_CSV"
fi

export MAP_ROOT REF_FA

python3 - <<'PY'
import os, csv
from pathlib import Path

MAP_ROOT   = Path(os.environ["MAP_ROOT"])
THRESHOLDS = [0.1, 1, 2, 5, 10, 20]

lorenz_header = ["Sample","BinLabel","BinRank","BinsTotal","BinFraction",
                 "CumulativeBinFraction","Depth","CumulativeDepth","CumulativeDepthFraction"]
for bl in ["10kb","50kb"]:
    f = MAP_ROOT / f"master_lorenz_curve_{bl}.csv"
    if not f.exists():
        with f.open("w", newline="") as fh:
            csv.writer(fh).writerow(lorenz_header)

bin_qc_header = [
    "Sample","BinLabel","BinsTotal","GC_valid_frac","DropoutFrac_depth0",
    "MeanDepth_all","MedianDepth_all","MeanDepth_nonzero","MedianDepth_nonzero",
    "IQR_nonzero","MAD_nonzero","CV_nonzero","Gini_nonzero",
    "Top1pct_share_nonzero","Top5pct_share_nonzero","Top10pct_share_nonzero",
    "MAPD_like_nonzero","GC_Spearman_nonzero","GC_Pearson_nonzero",
] + [f"Breadth_ge{t}".replace(".","p") for t in THRESHOLDS]
for bl in ["10kb","50kb"]:
    f = MAP_ROOT / f"master_bin_qc_metrics_{bl}.csv"
    if not f.exists():
        with f.open("w", newline="") as fh:
            csv.writer(fh).writerow(bin_qc_header)

picard_header = ["Sample",
    "READ_PAIRS_EXAMINED","UNPAIRED_READS_EXAMINED","UNMAPPED_READS",
    "READ_PAIR_DUPLICATES","UNPAIRED_READ_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE",
    "WINDOW_SIZE","TOTAL_CLUSTERS","AT_DROPOUT","GC_DROPOUT","MEAN_BIAS",
    "GENOME_TERRITORY","MEAN_COVERAGE","SD_COVERAGE","MEDIAN_COVERAGE","MAD_COVERAGE",
    "PCT_EXC_ADAPTER","PCT_EXC_MAPQ","PCT_EXC_DUPE","PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP","PCT_EXC_CAPPED","PCT_EXC_TOTAL",
    "PCT_0X","PCT_1X","PCT_5X","PCT_10X","PCT_20X","PCT_30X","PCT_50X","PCT_100X"]
f = MAP_ROOT / "master_picard_metrics.csv"
if not f.exists():
    with f.open("w", newline="") as fh:
        csv.writer(fh).writerow(picard_header)

print("Master CSV headers initialised.")
PY

# ============================================================
# LOOP 1: BWA mapping + samtools QC + markdup + metrics
# refmap.sorted.bam is KEPT — needed for normalization in LOOP 2.
# refmap.depth.txt and refmap.markdup.bam deleted immediately
# after metrics are extracted.
# ============================================================

shopt -s nullglob

for SAMPLE_DIR in "$SAMPLES_ROOT"/*_sc-*; do
  [ -d "$SAMPLE_DIR" ] || continue

  SAMPLE=$(basename "$SAMPLE_DIR")
  R1="$SAMPLE_DIR/${SAMPLE}_R1.fastq.gz"
  R2="$SAMPLE_DIR/${SAMPLE}_R2.fastq.gz"

  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "WARNING: missing R1/R2 for $SAMPLE — skipping."
    continue
  fi

  OUT_DIR="$MAP_ROOT/$SAMPLE"
  mkdir -p "$OUT_DIR"

  SCRATCH="$SCRATCH_BASE/$SAMPLE"
  mkdir -p "$SCRATCH"

  LOG_S1="$OUT_DIR/run.step1_mapping.log"

  echo "===== $(date) START LOOP1 $SAMPLE =====" | tee -a "$LOG_S1"

  if [ -f "$OUT_DIR/refmap.sorted.bam" ] && [ -f "$OUT_DIR/refmap.sorted.bam.bai" ]; then
    echo "INFO: $SAMPLE already has refmap.sorted.bam — skipping mapping." | tee -a "$LOG_S1"
  else
    if ! (
      SORT_TMP="$SCRATCH/${SAMPLE}.refmap.sorted.bam"
      TMP_PREFIX="$SCRATCH/${SAMPLE}.sorttmp"

      singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$BWA_SIF" \
        bwa mem -t "$THREADS" "$REF_FA" "$R1" "$R2" \
        2> "$OUT_DIR/bwa_mem.stderr.log" \
      | singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
          samtools view -b - \
      | singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
          samtools sort -@ "$THREADS" -T "$TMP_PREFIX" -o "$SORT_TMP" -

      cp -f "$SORT_TMP" "$OUT_DIR/refmap.sorted.bam"
      rm -f "$SORT_TMP"
      singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$OUT_DIR/refmap.sorted.bam"
    ) >> "$LOG_S1" 2>&1; then
      echo "ERROR: BWA mapping failed for $SAMPLE — see $LOG_S1. Continuing." | tee -a "$LOG_S1"
      continue
    fi
  fi

  if [ ! -f "$OUT_DIR/refmap.metrics_summary.csv" ]; then
    if ! (
      singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools flagstat \
        "$OUT_DIR/refmap.sorted.bam" > "$OUT_DIR/refmap.flagstat.txt"

      singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools depth -a \
        "$OUT_DIR/refmap.sorted.bam" > "$OUT_DIR/refmap.depth.txt"

      NAME_BAM="$SCRATCH/${SAMPLE}.name.bam"
      FIXMATE_BAM="$SCRATCH/${SAMPLE}.fixmate.bam"
      FIXSORT_BAM="$SCRATCH/${SAMPLE}.fixmate.sorted.bam"
      MD_BAM_TMP="$SCRATCH/${SAMPLE}.markdup.bam"
      TMP_PREFIX="$SCRATCH/${SAMPLE}.sorttmp"

      singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
        samtools sort -n -@ "$THREADS" -T "$TMP_PREFIX.name" -o "$NAME_BAM" "$OUT_DIR/refmap.sorted.bam"

      singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
        samtools fixmate -m "$NAME_BAM" "$FIXMATE_BAM"

      singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
        samtools sort -@ "$THREADS" -T "$TMP_PREFIX.fix" -o "$FIXSORT_BAM" "$FIXMATE_BAM"

      singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
        samtools markdup -@ "$THREADS" -s "$FIXSORT_BAM" "$MD_BAM_TMP" \
        2> "$OUT_DIR/markdup.metrics.txt"

      cp -f "$MD_BAM_TMP" "$OUT_DIR/refmap.markdup.bam"
      rm -f "$MD_BAM_TMP"
      singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$OUT_DIR/refmap.markdup.bam"
      singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools flagstat \
        "$OUT_DIR/refmap.markdup.bam" > "$OUT_DIR/refmap.markdup.flagstat.txt"

      FLAGSTAT="$OUT_DIR/refmap.flagstat.txt"
      DEPTH="$OUT_DIR/refmap.depth.txt"
      OUT_TSV="$OUT_DIR/refmap.metrics_summary.tsv"
      OUT_CSV="$OUT_DIR/refmap.metrics_summary.csv"

      TOTAL_READS=$(awk '/in total/ {print $1; exit}' "$FLAGSTAT")
      MAPPED_READS=$(awk '/ mapped \(/ && !/primary/ {print $1; exit}' "$FLAGSTAT")
      MAPPED_PCT=$(awk '/ mapped \(/ && !/primary/ {match($0,/\(([0-9.]+)%/,a); print (a[1]=="" ? "NA" : a[1]); exit}' "$FLAGSTAT")
      PAIRED=$(awk '/ paired in sequencing/ {print $1; exit}' "$FLAGSTAT")
      PROPER_PAIRED=$(awk '/ properly paired/ {print $1; exit}' "$FLAGSTAT")
      PROPER_PAIRED_PCT=$(awk '/ properly paired/ {match($0,/\(([0-9.]+)%/,a); print (a[1]=="" ? "NA" : a[1]); exit}' "$FLAGSTAT")
      SINGLETONS=$(awk '/ singletons/ {print $1; exit}' "$FLAGSTAT")
      MATE_DIFF_CHR=$(awk '/with mate mapped to a different chr$/ {print $1; exit}' "$FLAGSTAT")
      MATE_DIFF_CHR_MQ5=$(awk '/with mate mapped to a different chr \(mapQ>=5\)/ {print $1; exit}' "$FLAGSTAT")

      read -r POSITIONS MEAN_DEPTH_ALL MEAN_DEPTH_COV \
              COVER_GT0 FRAC_GT0 COVER_GE5 FRAC_GE5 \
              COVER_GE10 FRAC_GE10 ZERO_COUNT ZERO_FRAC <<EOF
$(awk '{d=$3;t++;sum+=d;if(d==0)z++;if(d>0){c0++;sumcov+=d}if(d>=5)c5++;if(d>=10)c10++;}
  END{printf("%d %.6f %.6f %d %.6f %d %.6f %d %.6f %d %.6f\n",
             t,sum/t,(c0?sumcov/c0:0),c0,c0/t,c5,c5/t,c10,c10/t,z,z/t)}' "$DEPTH")
EOF

      DUP_READS=$(awk '/ duplicates/ && !/primary/ {print $1; exit}' "$OUT_DIR/refmap.markdup.flagstat.txt" || true)

      {
        echo -e "Metric\tValue"
        echo -e "SAMPLE\t$SAMPLE"
        echo -e "MAP_DIR\t$OUT_DIR"
        echo -e "BAM\t$OUT_DIR/refmap.sorted.bam"
        echo -e "FLAGSTAT\t$FLAGSTAT"
        echo -e "DEPTH\t$DEPTH"
        echo -e "TOTAL_READS\t$TOTAL_READS"
        echo -e "MAPPED_READS\t$MAPPED_READS"
        echo -e "MAPPED_PCT\t$MAPPED_PCT"
        echo -e "PAIRED_IN_SEQUENCING\t$PAIRED"
        echo -e "PROPERLY_PAIRED\t$PROPER_PAIRED"
        echo -e "PROPERLY_PAIRED_PCT\t$PROPER_PAIRED_PCT"
        echo -e "SINGLETONS\t$SINGLETONS"
        echo -e "MATE_DIFF_CHR\t$MATE_DIFF_CHR"
        echo -e "MATE_DIFF_CHR_MAPQ_GE5\t$MATE_DIFF_CHR_MQ5"
        echo -e "POSITIONS_IN_DEPTH\t$POSITIONS"
        echo -e "MEAN_DEPTH_ALL_POS\t$MEAN_DEPTH_ALL"
        echo -e "MEAN_DEPTH_COVERED_POS\t$MEAN_DEPTH_COV"
        echo -e "FRAC_COVERED_GT0\t$FRAC_GT0"
        echo -e "FRAC_COVERED_GE5\t$FRAC_GE5"
        echo -e "FRAC_COVERED_GE10\t$FRAC_GE10"
        echo -e "FRAC_ZERO_DEPTH\t$ZERO_FRAC"
        echo -e "DUP_READS_MARKDUP_BAM\t$DUP_READS"
      } > "$OUT_TSV"

      awk -F'\t' 'BEGIN{OFS=","} NR==1{print "Metric","Value";next} {print $1,$2}' "$OUT_TSV" > "$OUT_CSV"

      echo "$SAMPLE,$MAPPED_PCT,$PROPER_PAIRED_PCT,$FRAC_GT0,$FRAC_GE5,$FRAC_GE10,$ZERO_FRAC,$MEAN_DEPTH_ALL,$MEAN_DEPTH_COV,$DUP_READS,$TOTAL_READS,$MAPPED_READS" \
        >> "$MASTER_REFMAP_CSV"
    ) >> "$LOG_S1" 2>&1; then
      echo "ERROR: Step 1B failed for $SAMPLE — see $LOG_S1. Continuing." | tee -a "$LOG_S1"
      continue
    fi
  fi

  rm -f "$OUT_DIR/refmap.depth.txt"
  rm -f "$OUT_DIR/refmap.markdup.bam" "$OUT_DIR/refmap.markdup.bam.bai"
  rm -f "$SCRATCH/${SAMPLE}.name.bam" "$SCRATCH/${SAMPLE}.fixmate.bam" \
        "$SCRATCH/${SAMPLE}.fixmate.sorted.bam"

  echo "===== $(date) DONE LOOP1 $SAMPLE =====" | tee -a "$LOG_S1"
done

echo "LOOP1 complete. Master table: $MASTER_REFMAP_CSV"

# ============================================================
# STEP 1.5: Compute normalization plan
# ============================================================

NORMALIZATION_PLAN="$MAP_ROOT/normalization_plan.csv"
echo "Sample,MAPPED_PCT,MAPPED_READS,KEEP_FOR_NORMALIZATION,NORMALIZATION_TARGET,NORMALIZATION_FRACTION,NORMALIZED_BAM" \
  > "$NORMALIZATION_PLAN"

TARGET_MAPPED_READS=$(awk -F',' -v minpct="$MIN_MAPPED_PCT" -v minreads="$MIN_MAPPED_READS" '
  NR>1 && $2!="NA" && $12!="NA" && ($2+0)>=minpct && ($12+0)>=minreads {print $12}
' "$REFERENCE_MASTER" | sort -n | head -1)

if [ -z "$TARGET_MAPPED_READS" ]; then
  echo "ERROR: No samples in REFERENCE_MASTER passed thresholds (MAPPED_PCT>=${MIN_MAPPED_PCT}%, MAPPED_READS>=${MIN_MAPPED_READS})"
  echo "ERROR: REFERENCE_MASTER=$REFERENCE_MASTER"
  exit 1
fi

echo "INFO: Normalization target from REFERENCE_MASTER = $TARGET_MAPPED_READS mapped reads"

awk -F',' -v OFS=',' \
    -v minpct="$MIN_MAPPED_PCT" \
    -v minreads="$MIN_MAPPED_READS" \
    -v target="$TARGET_MAPPED_READS" '
NR>1 {
  keep = ((($2+0) >= minpct) && (($12+0) >= minreads)) ? "YES" : "NO"
  frac = (keep=="YES" && ($12+0)>0) ? target/($12+0) : "NA"
  bam  = (keep=="YES") ? $1 "/refmap.normalized.sorted.bam" : "NA"
  print $1,$2,$12,keep,target,frac,bam
}' "$MASTER_REFMAP_CSV" >> "$NORMALIZATION_PLAN"

echo "INFO: Normalization plan written: $NORMALIZATION_PLAN"

# ============================================================
# LOOP 2: Normalize → Picard → mosdepth → per-sample metrics
# ============================================================

for SAMPLE_DIR in "$SAMPLES_ROOT"/*_sc-*; do
  [ -d "$SAMPLE_DIR" ] || continue

  SAMPLE=$(basename "$SAMPLE_DIR")
  OUT_DIR="$MAP_ROOT/$SAMPLE"

  [ -d "$OUT_DIR" ] || { echo "WARNING: $SAMPLE has no output dir — skipping LOOP2."; continue; }

  SCRATCH="$SCRATCH_BASE/$SAMPLE"
  mkdir -p "$SCRATCH"

  LOG_S15="$OUT_DIR/run.step1_5_normalization.log"
  LOG_S2="$OUT_DIR/run.step2_picard_mosdepth.log"
  LOG_S3="$OUT_DIR/run.step3_per_sample_metrics.log"

  NORM_BAM="$OUT_DIR/refmap.normalized.sorted.bam"
  RGBAM="$OUT_DIR/refmap.normalized.rg.bam"
  MDBAM="$OUT_DIR/refmap.normalized.picard_md.bam"

  echo "===== $(date) START STEP1.5 $SAMPLE =====" | tee -a "$LOG_S15"

  KEEP=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $4; exit}' "$NORMALIZATION_PLAN")
  FRACTION=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $6; exit}' "$NORMALIZATION_PLAN")

  if [ "$KEEP" != "YES" ]; then
    echo "INFO: $SAMPLE failed normalization thresholds — skipping LOOP2." | tee -a "$LOG_S15"
    rm -f "$OUT_DIR/refmap.sorted.bam" "$OUT_DIR/refmap.sorted.bam.bai"
    continue
  fi

  IN_BAM="$OUT_DIR/refmap.sorted.bam"

  if [ -f "$NORM_BAM" ] && [ -f "${NORM_BAM}.bai" ]; then
    echo "INFO: normalized BAM already exists for $SAMPLE — skipping normalization." | tee -a "$LOG_S15"
  else
    if [ ! -f "$IN_BAM" ]; then
      echo "WARNING: missing input BAM for normalization: $IN_BAM — skipping $SAMPLE." | tee -a "$LOG_S15"
      continue
    fi

    if ! (
      if awk -v x="$FRACTION" 'BEGIN{exit !(x > 1)}'; then
        echo "WARNING: $SAMPLE has fewer mapped reads than target; copying original BAM." | tee -a "$LOG_S15"
        cp -f "$IN_BAM" "$NORM_BAM"
        singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$NORM_BAM"
      elif awk -v x="$FRACTION" 'BEGIN{exit !(x >= 0.999999)}'; then
        cp -f "$IN_BAM" "$NORM_BAM"
        [ -f "${IN_BAM}.bai" ] && cp -f "${IN_BAM}.bai" "${NORM_BAM}.bai" || true
        singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$NORM_BAM"
      else
        TMP_NORM="$SCRATCH/${SAMPLE}.normalized.unsorted.bam"
        singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
          samtools view -@ "$THREADS" -s "${NORMALIZATION_SEED}.${FRACTION#0.}" -b "$IN_BAM" > "$TMP_NORM"
        singularity exec -B "$SCRATCH:$SCRATCH" -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" \
          samtools sort -@ "$THREADS" -T "$SCRATCH/${SAMPLE}.normsort" -o "$NORM_BAM" "$TMP_NORM"
        rm -f "$TMP_NORM"
        singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$NORM_BAM"
      fi
    ) >> "$LOG_S15" 2>&1; then
      echo "ERROR: normalization failed for $SAMPLE — see $LOG_S15. Continuing." | tee -a "$LOG_S15"
      continue
    fi
  fi

  rm -f "$IN_BAM" "${IN_BAM}.bai"

  echo "===== $(date) DONE STEP1.5 $SAMPLE =====" | tee -a "$LOG_S15"

  echo "===== $(date) START STEP2 $SAMPLE =====" | tee -a "$LOG_S2"

  if [ -f "$OUT_DIR/mosdepth_10kb.regions.bed.gz" ] && [ -f "$OUT_DIR/mosdepth_50kb.regions.bed.gz" ]; then
    echo "INFO: $SAMPLE already has mosdepth regions — skipping STEP2." | tee -a "$LOG_S2"
  else
    if ! (
      if [ ! -f "${NORM_BAM}.bai" ]; then
        singularity exec -B "$OUT_DIR:$OUT_DIR" "$SAMTOOLS_SIF" samtools index "$NORM_BAM"
      fi

      if [ ! -f "$RGBAM" ]; then
        java -jar "$PICARD_JAR" AddOrReplaceReadGroups \
          I="$NORM_BAM" O="$RGBAM" \
          RGID="$SAMPLE" RGLB="lib1" RGPL="ILLUMINA" \
          RGPU="${SAMPLE}.unit1" RGSM="$SAMPLE" \
          SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
      fi

      if [ ! -f "$MDBAM" ]; then
        java -jar "$PICARD_JAR" MarkDuplicates \
          I="$RGBAM" O="$MDBAM" M="$OUT_DIR/picard.MarkDuplicates.metrics.txt" \
          CREATE_INDEX=true ASSUME_SORT_ORDER=coordinate \
          VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null
      fi

      rm -f "$RGBAM" "${RGBAM%.bam}.bai" "${RGBAM}.bai"

      java -jar "$PICARD_JAR" CollectGcBiasMetrics \
        I="$MDBAM" O="$OUT_DIR/picard.gc_bias_metrics.txt" \
        CHART="$OUT_DIR/picard.gc_bias.pdf" \
        S="$OUT_DIR/picard.gc_bias_summary.txt" R="$REF_FA"

      java -jar "$PICARD_JAR" CollectWgsMetrics \
        I="$MDBAM" O="$OUT_DIR/picard.wgs_metrics.txt" R="$REF_FA"

      singularity exec "$MOSDEPTH_SIF" mosdepth \
        --threads "$THREADS" --by "$BIN_10KB" \
        "$OUT_DIR/mosdepth_10kb" "$MDBAM"

      singularity exec "$MOSDEPTH_SIF" mosdepth \
        --threads "$THREADS" --by "$BIN_50KB" \
        "$OUT_DIR/mosdepth_50kb" "$MDBAM"
    ) >> "$LOG_S2" 2>&1; then
      echo "ERROR: Step 2 failed for $SAMPLE — see $LOG_S2. Continuing." | tee -a "$LOG_S2"
      continue
    fi
  fi

  rm -f "$NORM_BAM" "${NORM_BAM}.bai"
  rm -f "$RGBAM" "${RGBAM%.bam}.bai" "${RGBAM}.bai"

  echo "===== $(date) DONE STEP2 $SAMPLE =====" | tee -a "$LOG_S2"

  echo "===== $(date) START STEP3 $SAMPLE =====" | tee -a "$LOG_S3"

  if [ -f "$OUT_DIR/bin_qc_10kb.tsv.gz" ] && [ -f "$OUT_DIR/bin_qc_50kb.tsv.gz" ] && \
     [ -f "$OUT_DIR/lorenz_curve_10kb.csv" ] && [ -f "$OUT_DIR/lorenz_curve_50kb.csv" ]; then
    echo "INFO: $SAMPLE already has bin_qc + lorenz outputs — skipping STEP3." | tee -a "$LOG_S3"
  else
    export OUT_DIR SAMPLE
    if ! python3 - >> "$LOG_S3" 2>&1 <<'PYEOF'
import os, gzip, math, csv, re
from pathlib import Path
from statistics import median

MAP_ROOT = Path(os.environ["MAP_ROOT"])
REF_FA   = Path(os.environ["REF_FA"])
OUT_DIR  = Path(os.environ["OUT_DIR"])
SAMPLE   = os.environ["SAMPLE"]

BIN_LABELS = ["10kb", "50kb"]
THRESHOLDS = [0.1, 1, 2, 5, 10, 20]

def load_fasta_as_dict(fa_path):
    seqs = {}; name = None; chunks = []
    with fa_path.open("rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                if name is not None: seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
        if name is not None: seqs[name] = "".join(chunks).upper()
    return seqs

def gc_fraction(seq):
    if not seq: return float("nan")
    seq = seq.upper()
    g=seq.count("G"); c=seq.count("C"); a=seq.count("A"); t=seq.count("T")
    d=a+t+g+c
    return (g+c)/d if d else float("nan")

def quantile_sorted(xs, q):
    n=len(xs)
    if n==0: return float("nan")
    if n==1: return xs[0]
    pos=(n-1)*q; lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    if lo==hi: return xs[lo]
    return xs[lo]*(1-(pos-lo))+xs[hi]*(pos-lo)

def mad(xs):
    if not xs: return float("nan")
    m=median(xs); return median([abs(x-m) for x in xs])

def gini(values):
    vals=[v for v in values if v is not None and not math.isnan(v)]
    if not vals: return float("nan")
    if all(v==0 for v in vals): return 0.0
    vals=sorted(vals); n=len(vals); s=sum(vals)
    if s==0: return 0.0
    cum=sum(i*x for i,x in enumerate(vals,1))
    return (2.0*cum)/(n*s)-(n+1.0)/n

def rankdata(a):
    n=len(a); order=sorted(range(n),key=lambda i:a[i])
    ranks=[0.0]*n; i=0; r=1
    while i<n:
        j=i
        while j+1<n and a[order[j+1]]==a[order[i]]: j+=1
        avg=(r+r+(j-i))/2.0
        for k in range(i,j+1): ranks[order[k]]=avg
        r+=(j-i+1); i=j+1
    return ranks

def pearson(x,y):
    n=len(x)
    if n<3: return float("nan")
    mx=sum(x)/n; my=sum(y)/n
    vx=sum((xi-mx)**2 for xi in x); vy=sum((yi-my)**2 for yi in y)
    if vx==0 or vy==0: return float("nan")
    return sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/math.sqrt(vx*vy)

def spearman(x,y): return pearson(rankdata(x),rankdata(y)) if len(x)>=3 else float("nan")

def top_share(vals,pct):
    if not vals: return float("nan")
    s=sum(vals)
    if s==0: return 0.0
    k=max(1,int(math.ceil(len(vals)*pct)))
    return sum(sorted(vals,reverse=True)[:k])/s

def mapd_like(depths_by_contig,eps=1e-3):
    diffs=[]
    for ds in depths_by_contig.values():
        prev=None
        for d in ds:
            if d<=0: prev=None; continue
            v=math.log(d+eps,2)
            if prev is not None: diffs.append(abs(v-prev))
            prev=v
    return median(diffs) if diffs else float("nan")

def parse_picard_first_table(path):
    if not path.exists(): return {}
    lines=[ln.rstrip("\n") for ln in path.open("rt",errors="replace")]
    i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line or line.startswith("#"): i+=1; continue
        header=re.split(r"\t+",line)
        if len(header)<2: header=re.split(r"\s+",line)
        if len(header)<2: i+=1; continue
        j=i+1
        while j<len(lines) and (not lines[j].strip() or lines[j].lstrip().startswith("#")): j+=1
        if j>=len(lines): return {}
        values=re.split(r"\t+",lines[j].strip())
        if len(values)!=len(header): values=re.split(r"\s+",lines[j].strip())
        if len(values)!=len(header): i+=1; continue
        return dict(zip(header,values))
    return {}

def fmt(v): return f"{v:.6f}" if isinstance(v,float) and not math.isnan(v) else "NA"

print(f"Loading reference FASTA: {REF_FA}")
ref = load_fasta_as_dict(REF_FA)
print(f"Loaded {len(ref)} contigs")

for bin_label in BIN_LABELS:
    regions_path = OUT_DIR / f"mosdepth_{bin_label}.regions.bed.gz"
    if not regions_path.exists():
        print(f"WARNING: missing {regions_path} — skipping {bin_label} for {SAMPLE}")
        continue

    lorenz_csv    = OUT_DIR / f"lorenz_curve_{bin_label}.csv"
    master_lorenz = MAP_ROOT / f"master_lorenz_curve_{bin_label}.csv"
    out_bins      = OUT_DIR / f"bin_qc_{bin_label}.tsv.gz"
    summary_tsv   = OUT_DIR / f"bin_qc_summary_{bin_label}.tsv"
    master_bq     = MAP_ROOT / f"master_bin_qc_metrics_{bin_label}.csv"

    bins_total=0; depth0=0; gc_valid=0
    ge_counts={t:0 for t in THRESHOLDS}
    depths_all=[]; depths_nz=[]; gc_nz=[]; depth_nz_for_gc=[]
    depths_by_contig={}; bin_rows=[]; lorenz_depths=[]

    with gzip.open(regions_path,"rt") as f_in:
        for line in f_in:
            parts=line.strip().split("\t")
            if len(parts)<4: continue
            contig,start_s,end_s,depth_s=parts[0],parts[1],parts[2],parts[3]
            start=int(start_s); end=int(end_s); d=float(depth_s)
            lorenz_depths.append(d); bins_total+=1; depths_all.append(d)
            depths_by_contig.setdefault(contig,[]).append(d)
            if d==0: depth0+=1
            for t in THRESHOLDS:
                if d>=t: ge_counts[t]+=1
            seq=ref.get(contig,"")
            gc=gc_fraction(seq[start:end]) if seq and end<=len(seq) else float("nan")
            if not math.isnan(gc): gc_valid+=1
            if d>0:
                depths_nz.append(d)
                if not math.isnan(gc): depth_nz_for_gc.append(d); gc_nz.append(gc)
            gc_str=f"{gc:.6f}" if not math.isnan(gc) else "NA"
            bin_rows.append(f"{contig}\t{start}\t{end}\t{end-start}\t{d:.6f}\t{gc_str}\n")

    if bins_total==0:
        print(f"WARNING: no bins in {regions_path} for {SAMPLE}"); continue

    lorenz_depths.sort()
    n=len(lorenz_depths); total_depth=sum(lorenz_depths)
    lorenz_hdr=["Sample","BinLabel","BinRank","BinsTotal","BinFraction",
                "CumulativeBinFraction","Depth","CumulativeDepth","CumulativeDepthFraction"]
    lorenz_rows=[]; cum_d=0.0
    for i,d in enumerate(lorenz_depths,1):
        cum_d+=d
        lorenz_rows.append([SAMPLE,bin_label,i,n,
                            f"{1.0/n:.10f}",f"{i/n:.10f}",
                            f"{d:.10f}",f"{cum_d:.10f}",
                            f"{(cum_d/total_depth if total_depth>0 else 0.0):.10f}"])
    with open(lorenz_csv,"w",newline="") as f:
        w=csv.writer(f); w.writerow(lorenz_hdr); w.writerows(lorenz_rows)
    with open(master_lorenz,"a",newline="") as f:
        csv.writer(f).writerows(lorenz_rows)
    print(f"Wrote Lorenz: {lorenz_csv}")

    with gzip.open(out_bins,"wt") as f_out:
        f_out.write("contig\tstart\tend\tlen\tmean_depth\tgc_frac\n")
        for row in bin_rows:
            f_out.write(row)
    print(f"Wrote bin_qc: {out_bins}")

    gc_valid_frac=gc_valid/bins_total; dropout_frac=depth0/bins_total
    mean_all=sum(depths_all)/bins_total; med_all=quantile_sorted(sorted(depths_all),0.5)
    breadths={t:ge_counts[t]/bins_total for t in THRESHOLDS}

    if depths_nz:
        mean_nz=sum(depths_nz)/len(depths_nz); nz_sorted=sorted(depths_nz)
        med_nz=quantile_sorted(nz_sorted,0.5)
        iqr_nz=quantile_sorted(nz_sorted,0.75)-quantile_sorted(nz_sorted,0.25)
        mad_nz=mad(depths_nz)
        cv_nz=(math.sqrt(sum((x-mean_nz)**2 for x in depths_nz)/(len(depths_nz)-1))/mean_nz
               if len(depths_nz)>=2 and mean_nz>0 else float("nan"))
        gini_nz=gini(depths_nz)
        top1=top_share(depths_nz,0.01); top5=top_share(depths_nz,0.05)
        top10=top_share(depths_nz,0.10); mapd=mapd_like(depths_by_contig)
    else:
        mean_nz=med_nz=iqr_nz=mad_nz=cv_nz=gini_nz=top1=top5=top10=mapd=float("nan")

    gc_spear=spearman(gc_nz,depth_nz_for_gc) if len(gc_nz)>=3 else float("nan")
    gc_pear=pearson(gc_nz,depth_nz_for_gc) if len(gc_nz)>=3 else float("nan")

    with summary_tsv.open("wt") as s:
        s.write("Metric\tValue\n")
        for k,v in [("Sample",SAMPLE),("BinLabel",bin_label),("BinsTotal",bins_total),
                    ("GC_valid_frac",f"{gc_valid_frac:.6f}"),("DropoutFrac_depth0",f"{dropout_frac:.6f}")]:
            s.write(f"{k}\t{v}\n")
        for t in THRESHOLDS:
            s.write(f"{('Breadth_ge'+str(t)).replace('.','p')}\t{breadths[t]:.6f}\n")
        for lbl,val in [("MeanDepth_all",fmt(mean_all)),("MedianDepth_all",fmt(med_all)),
                        ("MeanDepth_nonzero",fmt(mean_nz)),("MedianDepth_nonzero",fmt(med_nz)),
                        ("IQR_nonzero",fmt(iqr_nz)),("MAD_nonzero",fmt(mad_nz)),
                        ("CV_nonzero",fmt(cv_nz)),("Gini_nonzero",fmt(gini_nz)),
                        ("Top1pct_share_nonzero",fmt(top1)),("Top5pct_share_nonzero",fmt(top5)),
                        ("Top10pct_share_nonzero",fmt(top10)),("MAPD_like_nonzero",fmt(mapd)),
                        ("GC_Spearman_nonzero",fmt(gc_spear)),("GC_Pearson_nonzero",fmt(gc_pear))]:
            s.write(f"{lbl}\t{val}\n")

    bq_row=[SAMPLE,bin_label,bins_total,
        f"{gc_valid_frac:.6f}",f"{dropout_frac:.6f}",
        fmt(mean_all),fmt(med_all),fmt(mean_nz),fmt(med_nz),
        fmt(iqr_nz),fmt(mad_nz),fmt(cv_nz),fmt(gini_nz),
        fmt(top1),fmt(top5),fmt(top10),fmt(mapd),
        fmt(gc_spear),fmt(gc_pear)] + [f"{breadths[t]:.6f}" for t in THRESHOLDS]
    with open(master_bq,"a",newline="") as f:
        csv.writer(f).writerow(bq_row)

def pick(d,k): return d.get(k,"")
md  = parse_picard_first_table(OUT_DIR/"picard.MarkDuplicates.metrics.txt")
gc  = parse_picard_first_table(OUT_DIR/"picard.gc_bias_summary.txt")
wgs = parse_picard_first_table(OUT_DIR/"picard.wgs_metrics.txt")
picard_row=[SAMPLE,
    pick(md,"READ_PAIRS_EXAMINED"),pick(md,"UNPAIRED_READS_EXAMINED"),pick(md,"UNMAPPED_READS"),
    pick(md,"READ_PAIR_DUPLICATES"),pick(md,"UNPAIRED_READ_DUPLICATES"),
    pick(md,"PERCENT_DUPLICATION"),pick(md,"ESTIMATED_LIBRARY_SIZE"),
    pick(gc,"WINDOW_SIZE"),pick(gc,"TOTAL_CLUSTERS"),pick(gc,"AT_DROPOUT"),
    pick(gc,"GC_DROPOUT"),pick(gc,"MEAN_BIAS") or pick(gc,"MEAN_BIAS_COVERAGE"),
    pick(wgs,"GENOME_TERRITORY"),pick(wgs,"MEAN_COVERAGE"),pick(wgs,"SD_COVERAGE"),
    pick(wgs,"MEDIAN_COVERAGE"),pick(wgs,"MAD_COVERAGE"),
    pick(wgs,"PCT_EXC_ADAPTER"),pick(wgs,"PCT_EXC_MAPQ"),pick(wgs,"PCT_EXC_DUPE"),
    pick(wgs,"PCT_EXC_BASEQ"),pick(wgs,"PCT_EXC_OVERLAP"),pick(wgs,"PCT_EXC_CAPPED"),
    pick(wgs,"PCT_EXC_TOTAL"),pick(wgs,"PCT_0X"),pick(wgs,"PCT_1X"),pick(wgs,"PCT_5X"),
    pick(wgs,"PCT_10X"),pick(wgs,"PCT_20X"),pick(wgs,"PCT_30X"),
    pick(wgs,"PCT_50X"),pick(wgs,"PCT_100X")]
with open(MAP_ROOT/"master_picard_metrics.csv","a",newline="") as f:
    csv.writer(f).writerow(picard_row)
print(f"Appended Picard metrics for {SAMPLE}")
print("DONE per-sample metrics.")
PYEOF
    then
      echo "ERROR: per-sample Python metrics failed for $SAMPLE — see $LOG_S3. Continuing." | tee -a "$LOG_S3"
      continue
    fi
  fi

  rm -f "$MDBAM" "${MDBAM%.bam}.bai" "${MDBAM}.bai"
  rm -rf "$SCRATCH"

  echo "===== $(date) DONE ALL STEPS $SAMPLE =====" | tee -a "$LOG_S3"
done

echo "ALL STEPS COMPLETE."
echo "Refmap master CSV:          $MASTER_REFMAP_CSV"
echo "Normalization plan:         $NORMALIZATION_PLAN"
echo "Lorenz master CSVs:         $MAP_ROOT/master_lorenz_curve_10kb.csv"
echo "                            $MAP_ROOT/master_lorenz_curve_50kb.csv"
echo "Mosdepth bin QC masters:    $MAP_ROOT/master_bin_qc_metrics_10kb.csv"
echo "                            $MAP_ROOT/master_bin_qc_metrics_50kb.csv"
echo "Picard master CSV:          $MAP_ROOT/master_picard_metrics.csv"
