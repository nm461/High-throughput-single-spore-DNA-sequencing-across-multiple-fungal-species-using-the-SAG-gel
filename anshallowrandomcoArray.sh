#!/bin/bash
#$ -N ANrandCoAsm
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=1000:00:00
#$ -l s_vmem=16G,mem_req=16G
#$ -pe def_slot 8
#$ -t 1-100
#$ -o /home/nm461/normalisedassembly/250k/logs/ANshallowcoassemblyarray.out
#$ -e /home/nm461/normalisedassembly/250k/logs/ANshallowcoassemblyarray.log

set -euo pipefail
set -x

ITER_NUM="${SGE_TASK_ID:?SGE_TASK_ID not set}"

SUBSET_SIZE="${SUBSET_SIZE:-12}"
MIN_SAG="${MIN_SAG:-2}"
MAX_SAG="${MAX_SAG:-12}"
BASE_SEED="${BASE_SEED:-20250904}"
THREADS="${NSLOTS:-8}"

module use /usr/local/package/modulefiles
module load apptainer

SPADES_BIN="$HOME/bin/SPAdes-4.2.0-Linux/bin"
export PATH="$SPADES_BIN:$PATH"

command -v spades.py >/dev/null 2>&1 || {
  echo "ERROR: spades.py not found in PATH ($SPADES_BIN)"
  exit 1
}

OUT_ROOT="$HOME/normalisedassembly/250k/anshallow"
COA_DIR="$OUT_ROOT/ARRAYrandom_coassembly_bootstrap"
ITER_DIR="$COA_DIR/iterations"
META_DIR="$COA_DIR/metadata"

SIF_DIR="$HOME/sif_images"
BUSCO_LINEAGE="$HOME/fungi_odb10"

mkdir -p "$OUT_ROOT" "$COA_DIR" "$ITER_DIR" "$META_DIR"

SHARED_VENV="$OUT_ROOT/shared_venv"
if [ ! -d "$SHARED_VENV" ]; then
  python3 -m venv "$SHARED_VENV"
  source "$SHARED_VENV/bin/activate"
  pip install --upgrade pip
  pip install --no-cache-dir pandas numpy biopython
  deactivate
fi

export OUT_ROOT COA_DIR ITER_DIR META_DIR SIF_DIR BUSCO_LINEAGE
export SUBSET_SIZE MIN_SAG MAX_SAG BASE_SEED THREADS ITER_NUM

# ------------------------------------------------------------------
# Build eligible sample manifest once if missing
# ------------------------------------------------------------------
ELIGIBLE_CSV="$META_DIR/eligible_samples.csv"
if [ ! -s "$ELIGIBLE_CSV" ]; then
  source "$SHARED_VENV/bin/activate"
  python3 - << 'PY_ELIGIBLE'
import os
import csv
import pandas as pd

OUT_ROOT = os.environ["OUT_ROOT"]
META_DIR = os.environ["META_DIR"]
subset_size = int(os.environ["SUBSET_SIZE"])

master_csv = os.path.join(OUT_ROOT, "master_metrics.csv")
out_csv = os.path.join(META_DIR, "eligible_samples.csv")

if not os.path.exists(master_csv):
    raise SystemExit(f"ERROR: master_metrics.csv not found: {master_csv}")

df = pd.read_csv(master_csv)
if "Sample" not in df.columns:
    raise SystemExit("ERROR: 'Sample' column missing from master_metrics.csv")

samples = list(dict.fromkeys(df["Sample"].astype(str).tolist()))
rows = []

for sample in samples:
    sdir = os.path.join(OUT_ROOT, f"{sample}_out")
    if not os.path.isdir(sdir):
        continue
    if not os.path.exists(os.path.join(sdir, ".DONE")):
        continue

    candidates_r1 = [
        os.path.join(sdir, "norm_R1.fastq.gz"),
        os.path.join(sdir, "clean_R1.fastq.gz"),
        os.path.join(sdir, "trimmed_R1.fastq.gz"),
    ]
    candidates_r2 = [
        os.path.join(sdir, "norm_R2.fastq.gz"),
        os.path.join(sdir, "clean_R2.fastq.gz"),
        os.path.join(sdir, "trimmed_R2.fastq.gz"),
    ]

    r1 = next((p for p in candidates_r1 if os.path.exists(p)), None)
    r2 = next((p for p in candidates_r2 if os.path.exists(p)), None)

    if r1 and r2:
        rows.append((sample, r1, r2))

if len(rows) < subset_size:
    raise SystemExit(f"ERROR: only {len(rows)} eligible samples found, need {subset_size}")

tmp_csv = out_csv + ".tmp"
with open(tmp_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "r1", "r2"])
    w.writerows(rows)

os.replace(tmp_csv, out_csv)
print(f"Wrote eligible manifest: {out_csv}")
PY_ELIGIBLE
  deactivate
fi

ITER_NAME=$(printf "iter_%03d" "$ITER_NUM")
ITER_PATH="$ITER_DIR/$ITER_NAME"
SUBSET_TSV="$ITER_PATH/subset_manifest.tsv"
mkdir -p "$ITER_PATH"

export ITER_NAME ITER_PATH SUBSET_TSV

# ------------------------------------------------------------------
# Freeze random subset for this iteration
# ------------------------------------------------------------------
if [ ! -s "$SUBSET_TSV" ]; then
  source "$SHARED_VENV/bin/activate"
  python3 - << 'PY_SUBSET'
import os
import csv
import random
import pandas as pd

eligible_csv = os.path.join(os.environ["META_DIR"], "eligible_samples.csv")
iter_num = int(os.environ["ITER_NUM"])
subset_size = int(os.environ["SUBSET_SIZE"])
base_seed = int(os.environ["BASE_SEED"])
iter_path = os.environ["ITER_PATH"]

df = pd.read_csv(eligible_csv)
rng = random.Random(base_seed + iter_num)
picked = rng.sample(df.to_dict("records"), subset_size)

tmp = os.path.join(iter_path, "subset_manifest.tsv.tmp")
out = os.path.join(iter_path, "subset_manifest.tsv")

with open(tmp, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["order", "sample", "r1", "r2"])
    for i, row in enumerate(picked, start=1):
        w.writerow([i, row["sample"], row["r1"], row["r2"]])

os.replace(tmp, out)
print(f"Wrote {out}")
PY_SUBSET
  deactivate
fi

# ------------------------------------------------------------------
# Helper to parse BUSCO log into metrics CSV
# ------------------------------------------------------------------
write_busco_metrics_csv() {
  local sag="$1"
  local busco_log="$2"
  local out_csv="$3"

  export SAG_VALUE="$sag"
  export BUSCO_LOG="$busco_log"
  export BUSCO_OUT_CSV="$out_csv"

  source "$SHARED_VENV/bin/activate"
  python3 - << 'PY_PARSE'
import os, re, csv, math

sag = int(os.environ["SAG_VALUE"])
busco_log = os.environ["BUSCO_LOG"]
out_csv = os.environ["BUSCO_OUT_CSV"]

pat = re.compile(r"C:(\d+\.\d+)%.*D:(\d+\.\d+)%.*F:(\d+\.\d+)%.*M:(\d+\.\d+)%.*n:(\d+)")
comp = dup = frag = miss = float("nan")
n = float("nan")

with open(busco_log, encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        m = pat.search(line)
        if m:
            comp = float(m.group(1))
            dup = float(m.group(2))
            frag = float(m.group(3))
            miss = float(m.group(4))
            n = int(m.group(5))
            break

if math.isnan(comp):
    raise SystemExit(f"Could not parse BUSCO summary from {busco_log}")

tmp = out_csv + ".tmp"
with open(tmp, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["SAG", "Completeness", "Duplicated", "Fragmentation", "Missing", "BUSCO_n"])
    w.writerow([sag, comp, dup, frag, miss, n])

os.replace(tmp, out_csv)
print(f"Wrote {out_csv}")
PY_PARSE
  deactivate
}

# ------------------------------------------------------------------
# Run SAG 2..12 for this iteration
# Keep only: busco_metrics.csv, busco_<LABEL>.log, samples_used.tsv
# Everything else is deleted as soon as it is no longer needed.
# NOTE: if a SAG job is interrupted mid-run, delete its OUT_DIR to
#       restart cleanly, since intermediate files may be incomplete.
# ------------------------------------------------------------------
for K in $(seq "$MIN_SAG" "$MAX_SAG"); do
  LABEL=$(printf "SAG_%02d" "$K")
  OUT_DIR="$ITER_PATH/$LABEL"
  mkdir -p "$OUT_DIR"

  if [ -f "$OUT_DIR/.DONE" ]; then
    echo "[$ITER_NAME $LABEL] already done; skipping"
    continue
  fi

  export K_VALUE="$K"
  export OUT_DIR SUBSET_TSV ITER_PATH

  # Write per-SAG sample list and R1/R2 file lists
  source "$SHARED_VENV/bin/activate"
  python3 - << 'PY_LISTS'
import os, csv

subset_tsv = os.environ["SUBSET_TSV"]
out_dir = os.environ["OUT_DIR"]
k = int(os.environ["K_VALUE"])

with open(subset_tsv, encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

rows = rows[:k]
if len(rows) < k:
    raise SystemExit(f"Need {k} rows, found {len(rows)}")

with open(os.path.join(out_dir, "samples_used.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["order", "sample", "r1", "r2"])
    for row in rows:
        w.writerow([row["order"], row["sample"], row["r1"], row["r2"]])

with open(os.path.join(out_dir, "R1.list"), "w", encoding="utf-8") as f:
    for row in rows:
        f.write(row["r1"] + "\n")

with open(os.path.join(out_dir, "R2.list"), "w", encoding="utf-8") as f:
    for row in rows:
        f.write(row["r2"] + "\n")
PY_LISTS
  deactivate

  # --- Concatenate reads ---
  if [ ! -s "$OUT_DIR/R1_combined.fastq.gz" ] || [ ! -s "$OUT_DIR/R2_combined.fastq.gz" ]; then
    while read -r f; do zcat "$f"; done < "$OUT_DIR/R1.list" | gzip > "$OUT_DIR/R1_combined.fastq.gz"
    while read -r f; do zcat "$f"; done < "$OUT_DIR/R2.list" | gzip > "$OUT_DIR/R2_combined.fastq.gz"
  fi

  # --- SPAdes assembly ---
  if [ ! -s "$OUT_DIR/spades_out/contigs.fasta" ]; then
    spades.py \
      --sc \
      -1 "$OUT_DIR/R1_combined.fastq.gz" \
      -2 "$OUT_DIR/R2_combined.fastq.gz" \
      -o "$OUT_DIR/spades_out" \
      --threads "$THREADS" \
      --memory 64 \
      > "$OUT_DIR/spades.log" 2>&1
  fi

  ASSEMBLY="$OUT_DIR/spades_out/contigs.fasta"
  [ -s "$ASSEMBLY" ] || { echo "ERROR: SPAdes contigs missing: $ASSEMBLY"; exit 1; }

  # --- BWA map back to assembly (needed for Pilon polishing) ---
  if [ ! -s "$OUT_DIR/mapped.sorted.bam" ]; then
    singularity exec "$SIF_DIR/bwa.sif" bwa index "$ASSEMBLY"
    singularity exec "$SIF_DIR/bwa.sif" bwa mem -t "$THREADS" \
      "$ASSEMBLY" "$OUT_DIR/R1_combined.fastq.gz" "$OUT_DIR/R2_combined.fastq.gz" \
      | singularity exec "$SIF_DIR/samtools.sif" samtools view -bS - \
      > "$OUT_DIR/mapped.bam"
    singularity exec "$SIF_DIR/samtools.sif" samtools sort -@ "$THREADS" \
      -o "$OUT_DIR/mapped.sorted.bam" "$OUT_DIR/mapped.bam"
    singularity exec "$SIF_DIR/samtools.sif" samtools index "$OUT_DIR/mapped.sorted.bam"
    rm -f "$OUT_DIR/mapped.bam"  # unsorted no longer needed
  fi

  # Combined reads no longer needed (SPAdes + BWA both done)
  rm -f "$OUT_DIR/R1_combined.fastq.gz" "$OUT_DIR/R2_combined.fastq.gz"

  # --- Pilon polishing ---
  if [ ! -s "$OUT_DIR/pilon_corrected.fasta" ]; then
    singularity exec "$SIF_DIR/pilon.sif" java -Xmx60G -jar /pilon/pilon.jar \
      --genome "$ASSEMBLY" \
      --frags "$OUT_DIR/mapped.sorted.bam" \
      --output "$OUT_DIR/pilon_corrected" \
      --threads "$THREADS" \
      > "$OUT_DIR/pilon.log" 2>&1
  fi

  # BAM, BWA index, and SPAdes output no longer needed
  rm -f "$OUT_DIR/mapped.sorted.bam" "$OUT_DIR/mapped.sorted.bam.bai"
  rm -f "${ASSEMBLY}.amb" "${ASSEMBLY}.ann" "${ASSEMBLY}.bwt" "${ASSEMBLY}.pac" "${ASSEMBLY}.sa"
  rm -rf "$OUT_DIR/spades_out"

  # --- Seqkit: filter to contigs >= 1 kb ---
  FILTERED_ASSEMBLY="$OUT_DIR/pilon_corrected.fasta"
  if [ -f "$SIF_DIR/seqkit.sif" ]; then
    if [ ! -s "$OUT_DIR/pilon_corrected.1kb.fasta" ]; then
      singularity exec "$SIF_DIR/seqkit.sif" seqkit seq -m 1000 "$OUT_DIR/pilon_corrected.fasta" \
        > "$OUT_DIR/pilon_corrected.1kb.fasta"
    fi
    rm -f "$OUT_DIR/pilon_corrected.fasta"  # unfiltered no longer needed
    FILTERED_ASSEMBLY="$OUT_DIR/pilon_corrected.1kb.fasta"
  fi

  # --- BUSCO ---
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  BUSCO_LOG="$OUT_DIR/busco_${LABEL}.log"
  if [ ! -s "$BUSCO_LOG" ]; then
    singularity exec -B "$OUT_DIR:$OUT_DIR" "$SIF_DIR/busco_v6.0.0_cv1.sif" \
      busco \
      -i "$FILTERED_ASSEMBLY" \
      -o "busco_${LABEL}" \
      -l "$BUSCO_LINEAGE" \
      -m genome \
      --cpu 2 \
      --out_path "$OUT_DIR" \
      > "$BUSCO_LOG" 2>&1
  fi

  # Parse BUSCO log → metrics CSV (the two files we keep)
  write_busco_metrics_csv "$K" "$BUSCO_LOG" "$OUT_DIR/busco_metrics.csv"

  # Delete filtered assembly and BUSCO output directory; keep only log + metrics CSV
  rm -f "$FILTERED_ASSEMBLY"
  rm -rf "$OUT_DIR/busco_${LABEL}"

  touch "$OUT_DIR/.DONE"
  echo "[$ITER_NAME $LABEL] complete"
done

touch "$ITER_PATH/.DONE"
echo "[$ITER_NAME] complete"
