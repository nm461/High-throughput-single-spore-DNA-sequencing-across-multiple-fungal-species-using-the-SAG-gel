#!/bin/bash
#$ -N ANdeepNorm250k
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=400:00:00
#$ -l s_vmem=16G,mem_req=16G
#$ -pe def_slot 8
#$ -o logs/AndeepNorm250k.out
#$ -e logs/AndeepNorm250k.log

set -e
set -x

START_AT="${START_AT:-all}"   # all or coassembly

# =========================
# Modules
# =========================
module use /usr/local/package/modulefiles
module load apptainer
module load blast+/2.15.0
module load repeatmasker/4.1.6 || true   # harmless; RepeatMasker is skipped

export PATH="$HOME/tools/rmblast/rmblast-2.14.1/bin:$PATH"

# =========================
# Tools on PATH (SPAdes)
# =========================
SPADES_BIN="${SPADES_BIN:-$HOME/bin/SPAdes-4.2.0-Linux/bin}"
export PATH="$SPADES_BIN:$PATH"

# Fail early if SPAdes missing
command -v spades.py >/dev/null 2>&1 || { echo "ERROR: spades.py not found in PATH ($SPADES_BIN)"; exit 1; }

# =========================
# Directories & inputs
# =========================
OUT_ROOT="${OUT_ROOT:-$PWD/results/andeep_250k}"
SIF_DIR="${SIF_DIR:-$HOME/sif_images}"
DATA_DIR="${DATA_DIR:-$PWD/data/andeep_fastq}"
DB_DIR="${DB_DIR:-$HOME/databases/kraken2_plusPFP_16G/k2_pluspfp_16G}"   # Kraken2 DB dir (contains *.k2d)
BUSCO_LINEAGE="${BUSCO_LINEAGE:-$HOME/fungi_odb10}"                      # BUSCO lineage dir
HOST_REF="${HOST_REF:-$HOME/databases/human_masked/human_masked.fasta}"
MULTIQC_CUSTOM="${MULTIQC_CUSTOM:-$HOME/multiqc_custom}"
SM_DB_FUNGI="${SM_DB_FUNGI:-$HOME/databases/sourmash/fungi/genbank-2022.03-fungi-k21.zip}"
SM_DB_MAMMAL="${SM_DB_MAMMAL:-$HOME/databases/sourmash/mammals/hg38.sig.zip}"

# Read normalization (post-QC / post-host-decontam)
RASUSA_SIF="$SIF_DIR/rasusa.sif"
NORMALIZE_READS="${NORMALIZE_READS:-250000}"
MIN_READS="${MIN_READS:-250000}"

# Limit how many samples to run this time
MAX_SAMPLES="${MAX_SAMPLES:-16}"

mkdir -p "$OUT_ROOT"
export SIF_DIR
export OUT_ROOT

# =========================
# Shared Python venv (once)
# =========================
SHARED_VENV="$OUT_ROOT/shared_venv"
if [ ! -d "$SHARED_VENV" ]; then
  python3 -m venv "$SHARED_VENV"
  source "$SHARED_VENV/bin/activate"
  pip install --upgrade pip
  pip install --no-cache-dir pandas numpy matplotlib seaborn biopython scipy PyPDF2
  deactivate
fi

if [ "$START_AT" != "coassembly" ]; then

# =========================
# MAIN LOOP OVER SAMPLES (capped at MAX_SAMPLES)
# =========================
processed=0
shopt -s nullglob

for R1 in "$DATA_DIR"/AN*/*_R1.fastq.gz; do
  if [ "$processed" -ge "$MAX_SAMPLES" ]; then
    echo "Reached MAX_SAMPLES=$MAX_SAMPLES; stopping."
    break
  fi

  sample=$(basename "$R1" _R1.fastq.gz)
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

  if [ ! -f "$R2" ]; then
    echo "WARNING: missing R2 for sample $sample — skipping."
    continue
  fi

  processed=$((processed + 1))

  OUT_DIR="$OUT_ROOT/${sample}_out"
  mkdir -p "$OUT_DIR"

  # Skip if this sample already finished successfully
  if [ -f "$OUT_DIR/.DONE" ]; then
    echo "==== $sample already done; skipping ===="
    continue
  fi
  echo "==== Processing $sample ===="

  # ---------- 1. Fastp ----------
  singularity exec "$SIF_DIR/fastp.sif" fastp \
    -i "$R1" -I "$R2" \
    -o "$OUT_DIR/trimmed_R1.fastq.gz" \
    -O "$OUT_DIR/trimmed_R2.fastq.gz" \
    --thread "$NSLOTS" \
    --length_required 100 \
    --json "$OUT_DIR/fastp.json" \
    > "$OUT_DIR/fastp.log" 2>&1

  # ---------- 1.5 Optional: Host decontamination (BWA) ----------
  if [ -f "$HOST_REF" ]; then
      echo "Running host decontamination using BWA..."

      singularity exec "$SIF_DIR/bwa.sif" bwa mem -t "$NSLOTS" \
          "$HOST_REF" "$OUT_DIR/trimmed_R1.fastq.gz" "$OUT_DIR/trimmed_R2.fastq.gz" \
        | singularity exec "$SIF_DIR/samtools.sif" samtools view -bS - \
        > "$OUT_DIR/host_mapped.bam"

      singularity exec "$SIF_DIR/samtools.sif" samtools view -b -f 12 -F 256 \
          "$OUT_DIR/host_mapped.bam" > "$OUT_DIR/clean_unmapped.bam"

      singularity exec "$SIF_DIR/samtools.sif" samtools fastq \
          -1 "$OUT_DIR/clean_R1.fastq.gz" \
          -2 "$OUT_DIR/clean_R2.fastq.gz" \
          -0 /dev/null -s /dev/null -n \
          "$OUT_DIR/clean_unmapped.bam"

      R1_CLEAN="$OUT_DIR/clean_R1.fastq.gz"
      R2_CLEAN="$OUT_DIR/clean_R2.fastq.gz"
  else
      echo "Host reference not found — skipping host decontam."
      R1_CLEAN="$OUT_DIR/trimmed_R1.fastq.gz"
      R2_CLEAN="$OUT_DIR/trimmed_R2.fastq.gz"
  fi

  # ---------- 1.75 Read normalization (rasusa to 1,000,000 read pairs) ----------
  [ -f "$RASUSA_SIF" ] || { echo "ERROR: rasusa SIF not found at $RASUSA_SIF"; exit 1; }

  READS_AFTER_CLEAN=$(zcat "$R1_CLEAN" | awk 'END{print NR/4}')
  echo "[$sample] usable read pairs after QC/host removal: $READS_AFTER_CLEAN" | tee "$OUT_DIR/read_normalization.log"

  if [ "$READS_AFTER_CLEAN" -lt "$MIN_READS" ]; then
    echo "[$sample] SKIP: usable reads ($READS_AFTER_CLEAN) < MIN_READS ($MIN_READS)" | tee -a "$OUT_DIR/read_normalization.log"
    continue
  fi

  if [ "$READS_AFTER_CLEAN" -gt "$NORMALIZE_READS" ]; then
    echo "[$sample] Downsampling to $NORMALIZE_READS read pairs with rasusa" | tee -a "$OUT_DIR/read_normalization.log"

    singularity exec "$RASUSA_SIF" rasusa reads \
      -n "$NORMALIZE_READS" \
      -o "$OUT_DIR/norm_R1.fastq.gz" \
      -o "$OUT_DIR/norm_R2.fastq.gz" \
      "$R1_CLEAN" "$R2_CLEAN" \
      >> "$OUT_DIR/read_normalization.log" 2>&1

    R1_NORM="$OUT_DIR/norm_R1.fastq.gz"
    R2_NORM="$OUT_DIR/norm_R2.fastq.gz"
  else
    echo "[$sample] No downsampling needed; using cleaned reads directly" | tee -a "$OUT_DIR/read_normalization.log"
    R1_NORM="$R1_CLEAN"
    R2_NORM="$R2_CLEAN"
  fi

  # ---------- 2. SPAdes (single-cell mode) ----------
  export PATH="$SPADES_BIN:$PATH"
  spades.py \
    --sc \
    -1 "$R1_NORM" \
    -2 "$R2_NORM" \
    -o "$OUT_DIR/spades_out" \
    --threads "$NSLOTS" \
    --memory 120 \
    > "$OUT_DIR/spades.log" 2>&1

  ASSEMBLY="$OUT_DIR/spades_out/contigs.fasta"

  # ---------- 3. Kraken2 ----------
  KRAKEN_DIR="$OUT_DIR/contamination_checks/kraken"
  mkdir -p "$KRAKEN_DIR"

  singularity exec "$SIF_DIR/kraken2.sif" kraken2 \
    --db "$DB_DIR" \
    --report "$KRAKEN_DIR/kraken_report.txt" \
    --paired "$R1_NORM" "$R2_NORM" \
    --threads "$NSLOTS" \
    > "$KRAKEN_DIR/kraken_output.txt" 2>&1

  # ---------- 4. BWA/Samtools ----------
  singularity exec "$SIF_DIR/bwa.sif" bwa index "$ASSEMBLY"
  singularity exec "$SIF_DIR/bwa.sif" bwa mem -t "$NSLOTS" \
    "$ASSEMBLY" "$R1_NORM" "$R2_NORM" \
    | singularity exec "$SIF_DIR/samtools.sif" samtools view -bS - \
    > "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools sort -@ "$NSLOTS" \
    -o "$OUT_DIR/mapped.sorted.bam" "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools index "$OUT_DIR/mapped.sorted.bam"

  # ---------- 5. Pilon ----------
  singularity exec "$SIF_DIR/pilon.sif" java -Xmx120G -jar /pilon/pilon.jar \
    --genome "$ASSEMBLY" \
    --frags "$OUT_DIR/mapped.sorted.bam" \
    --output "$OUT_DIR/pilon_corrected" \
    --threads "$NSLOTS" \
    > "$OUT_DIR/pilon.log" 2>&1

  PILON_ASSEMBLY="$OUT_DIR/pilon_corrected.fasta"

  # ---------- 5.5 Contig filtering (>= 1000 bp, using seqkit SIF if available) ----------
  FILTERED_ASSEMBLY="$PILON_ASSEMBLY"
  SEQKIT_SIF="$SIF_DIR/seqkit.sif"

  if [ -f "$SEQKIT_SIF" ]; then
    echo "Filtering contigs < 1000 bp using seqkit (via SIF)..."
    singularity exec "$SEQKIT_SIF" seqkit seq -m 1000 "$PILON_ASSEMBLY" > "$OUT_DIR/pilon_corrected.1kb.fasta"
    FILTERED_ASSEMBLY="$OUT_DIR/pilon_corrected.1kb.fasta"
  else
    echo "Seqkit SIF not found at $SEQKIT_SIF — skipping contig length filter."
  fi

  # ---------- 6. BUSCO (on polished, filtered assembly, with diagnostics) ----------
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  TS=$(date +%Y%m%d%H%M%S)
  BUSCO_RUN="busco_pilon_${TS}"
  BUSCO_LOG="$OUT_DIR/${BUSCO_RUN}.log"

  {
    echo "=== BUSCO preflight for $sample ==="
    date
    echo "ASSEMBLY (Pilon-corrected, filtered if available): $FILTERED_ASSEMBLY"; ls -lh "$FILTERED_ASSEMBLY" || true
    echo "BUSCO_LINEAGE: $BUSCO_LINEAGE"; ls -lh "$BUSCO_LINEAGE" | head || true
    singularity exec -B "$OUT_DIR:$OUT_DIR" "$SIF_DIR/busco_v6.0.0_cv1.sif" bash -lc '
      echo "busco version:"; busco --version || true
      echo "pwd: $(pwd)"
    ' || true
  } > "$OUT_DIR/busco_diag.log" 2>&1

  singularity exec -B "$OUT_DIR:$OUT_DIR" "$SIF_DIR/busco_v6.0.0_cv1.sif" \
    busco \
      -i "$FILTERED_ASSEMBLY" \
      -o "$BUSCO_RUN" \
      -l "$BUSCO_LINEAGE" \
      -m genome \
      --cpu 2 \
      --out_path "$OUT_DIR" \
    > "$BUSCO_LOG" 2>&1

  # ---------- 7. Sourmash ----------
  SM_OUT="$OUT_DIR/contamination_checks/sourmash"
  mkdir -p "$SM_OUT"
  SIF_SM="$SIF_DIR/sourmash.sif"

  singularity exec "$SIF_SM" sourmash sketch dna -p k=21,scaled=1000 \
    -o "$SM_OUT/assembly.sig" "$FILTERED_ASSEMBLY"

  [ -f "$SM_DB_FUNGI" ]  || { echo "Missing SM_DB_FUNGI: $SM_DB_FUNGI"  >&2; exit 1; }
  [ -f "$SM_DB_MAMMAL" ] || { echo "Missing SM_DB_MAMMAL: $SM_DB_MAMMAL" >&2; exit 1; }

  singularity exec "$SIF_SM" sourmash gather \
    "$SM_OUT/assembly.sig" \
    "$SM_DB_FUNGI" \
    "$SM_DB_MAMMAL" \
    -o "$SM_OUT/gather_results.csv"

  # ---------- 8. Per-sample text report (PDF) + CSV metrics + GraphData CSVs ----------
  source "$SHARED_VENV/bin/activate"
  export OUT_DIR OUT_ROOT SIF_DIR sample
  python3 - << 'PY_PER_SAMPLE'
import os, json, re, glob, datetime, csv, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
ddf = None
kdf = None

OUT_DIR  = os.environ.get("OUT_DIR", "")
OUT_ROOT = os.environ.get("OUT_ROOT", "")
SIF_DIR  = os.path.expanduser(os.environ.get("SIF_DIR", "~/sif_images"))
SAMPLE   = os.environ.get("sample", "")

if not OUT_DIR or not OUT_ROOT or not SAMPLE:
    raise SystemExit("OUT_DIR / OUT_ROOT / sample not set")

def pick_font():
    preferred = [
        "DejaVu Sans", "Arial", "Liberation Sans",
        "Noto Sans", "Noto Sans CJK JP"
    ]
    paths = fm.findSystemFonts(fontext="ttf") + fm.findSystemFonts(fontext="ttc")
    names = {}
    for p in paths:
        try:
            names[fm.get_font(p).family_name] = p
        except Exception:
            pass
    for name in preferred:
        if name in names:
            return fm.FontProperties(fname=names[name])
    return None

FONT = pick_font()

def ensure_depth():
    depth = os.path.join(OUT_DIR, "depth.txt")
    bam   = os.path.join(OUT_DIR, "mapped.sorted.bam")
    sif   = os.path.join(SIF_DIR, "samtools.sif")
    if os.path.exists(depth):
        return depth
    if os.path.exists(bam) and os.path.exists(sif):
        with open(depth, "w") as out:
            subprocess.call(
                ["singularity", "exec", sif, "samtools", "depth", "-a", bam],
                stdout=out
            )
        return depth
    return None

fjson = os.path.join(OUT_DIR, "fastp.json")
q30 = None
r1len = None
r2len = None
pairs_before = None
pairs_after  = None

if os.path.exists(fjson):
    try:
        jj = json.load(open(fjson))
        summary = jj.get("summary", {})
        before  = summary.get("before_filtering", {})
        after   = summary.get("after_filtering", {})
        q30     = after.get("q30_rate")
        r1len   = after.get("read1_mean_length")
        r2len   = after.get("read2_mean_length")
        pairs_before = before.get("total_reads")
        pairs_after  = after.get("total_reads")
    except Exception:
        pass

fa_candidates = [
    os.path.join(OUT_DIR, "pilon_corrected.1kb.fasta"),
    os.path.join(OUT_DIR, "pilon_corrected.fasta"),
    os.path.join(OUT_DIR, "spades_out", "contigs.fasta"),
]
assembly_path = next((p for p in fa_candidates if os.path.exists(p)), None)

total_len = np.nan
n50 = np.nan
max_ctg = np.nan
n_contigs = np.nan

if assembly_path:
    lens = [len(rec.seq) for rec in SeqIO.parse(assembly_path, "fasta")]
    if lens:
        n_contigs = len(lens)
        total_len = float(sum(lens))
        max_ctg   = float(max(lens))
        s = sorted(lens, reverse=True)
        csum = 0
        half = total_len / 2.0
        for L in s:
            csum += L
            if csum >= half:
                n50 = float(L)
                break

busco_c = busco_d = busco_f = busco_m = np.nan
busco_n = np.nan

log_cands = sorted(
    glob.glob(os.path.join(OUT_DIR, "busco_pilon_*.log")),
    key=os.path.getmtime,
    reverse=True
)
if log_cands:
    blog = log_cands[0]
    for line in open(blog, encoding="utf-8", errors="ignore"):
        m = re.search(
            r"C:(\d+\.\d+)%.*D:(\d+\.\d+)%.*F:(\d+\.\d+)%.*M:(\d+\.\d+)%.*n:(\d+)",
            line
        )
        if m:
            c, d, f, m_, n = m.groups()
            busco_c = float(c)
            busco_d = float(d)
            busco_f = float(f)
            busco_m = float(m_)
            busco_n = int(n)
            break

depth_path = ensure_depth()
depth_mean = depth_median = np.nan
if depth_path and os.path.exists(depth_path):
    try:
        ddf = pd.read_csv(depth_path, sep="\t", header=None,
                          names=["contig", "pos", "depth"])
        depth_mean   = float(ddf["depth"].mean())
        depth_median = float(ddf["depth"].median())
    except Exception:
        pass

krep = os.path.join(OUT_DIR, "contamination_checks", "kraken", "kraken_report.txt")
kraken_top = []
if os.path.exists(krep):
    try:
        kdf = pd.read_csv(
            krep, sep="\t", header=None,
            names=["perc","reads","tax_reads","rank","taxid","name"]
        )
        kdf["name"] = kdf["name"].astype(str).str.strip()

        sdf = kdf[kdf["rank"] == "S"].copy()
        sdf = sdf.sort_values("perc", ascending=False)

        top = sdf.head(5)
        for _, r in top.iterrows():
            kraken_top.append((r["name"], float(r["perc"])))

        kraken_long_path = os.path.join(
            OUT_DIR, "contamination_checks", "kraken", "kraken_species_long.csv"
        )
        os.makedirs(os.path.dirname(kraken_long_path), exist_ok=True)
        out = sdf.head(50)[["name","perc","reads","tax_reads","taxid"]].copy()
        out.insert(0, "Sample", SAMPLE)
        out.to_csv(kraken_long_path, index=False)

    except Exception:
        kdf = None

sm_gather = os.path.join(
    OUT_DIR, "contamination_checks", "sourmash", "gather_results.csv"
)
sm_long_path = os.path.join(
    OUT_DIR, "contamination_checks", "sourmash", "sourmash_gather_long.csv"
)

def pick_sourmash_metric(df):
    for c in ["f_match", "f_unique_to_query", "f_unique_weighted", "f_query_match"]:
        if c in df.columns:
            return c
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            return c
    return None

if os.path.exists(sm_gather):
    try:
        sdf = pd.read_csv(sm_gather)
        metric = pick_sourmash_metric(sdf)
        name_col = (
            "name" if "name" in sdf.columns
            else ("match_name" if "match_name" in sdf.columns else None)
        )

        if metric and name_col:
            out = sdf[[name_col, metric]].copy()
            out.rename(columns={name_col: "match", metric: "value"}, inplace=True)
            out.insert(0, "Sample", SAMPLE)
            out = out.sort_values("value", ascending=False).head(50)
            os.makedirs(os.path.dirname(sm_long_path), exist_ok=True)
            out.to_csv(sm_long_path, index=False)
    except Exception:
        pass

if kraken_top:
    kt_path = os.path.join(OUT_DIR, "kraken_top5.csv")
    os.makedirs(os.path.dirname(kt_path), exist_ok=True)
    with open(kt_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Taxon", "Percent_reads"])
        for name, perc in kraken_top:
            w.writerow([name, perc])

sample_metrics_path = os.path.join(OUT_DIR, "sample_metrics.csv")
cols = [
    "Sample",
    "Assembly_path",
    "TotalLength_bp",
    "N50_bp",
    "MaxContig_bp",
    "NumContigs",
    "BUSCO_Complete",
    "BUSCO_Duplicated",
    "BUSCO_Fragmented",
    "BUSCO_Missing",
    "BUSCO_n",
    "Depth_mean",
    "Depth_median",
    "Fastp_Q30_rate",
    "Fastp_R1_mean_len",
    "Fastp_R2_mean_len",
    "Reads_before",
    "Reads_after",
]
row = {
    "Sample": SAMPLE,
    "Assembly_path": assembly_path or "",
    "TotalLength_bp": total_len,
    "N50_bp": n50,
    "MaxContig_bp": max_ctg,
    "NumContigs": n_contigs,
    "BUSCO_Complete": busco_c,
    "BUSCO_Duplicated": busco_d,
    "BUSCO_Fragmented": busco_f,
    "BUSCO_Missing": busco_m,
    "BUSCO_n": busco_n,
    "Depth_mean": depth_mean,
    "Depth_median": depth_median,
    "Fastp_Q30_rate": float(q30) if q30 is not None else np.nan,
    "Fastp_R1_mean_len": float(r1len) if r1len is not None else np.nan,
    "Fastp_R2_mean_len": float(r2len) if r2len is not None else np.nan,
    "Reads_before": int(pairs_before) if pairs_before is not None else np.nan,
    "Reads_after": int(pairs_after) if pairs_after is not None else np.nan,
}

with open(sample_metrics_path, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerow(row)

master_path = os.path.join(OUT_ROOT, "master_metrics.csv")
write_header = not os.path.exists(master_path)
with open(master_path, "a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    if write_header:
        w.writeheader()
    w.writerow(row)

GDIR = os.path.join(OUT_DIR, "GraphData")
os.makedirs(GDIR, exist_ok=True)

dirs = {
    "ReadQuality": os.path.join(GDIR, "ReadQuality"),
    "ContigLengths": os.path.join(GDIR, "ContigLengths"),
    "AssemblyStats": os.path.join(GDIR, "AssemblyStats"),
    "BUSCO": os.path.join(GDIR, "BUSCO"),
    "Coverage": os.path.join(GDIR, "Coverage"),
    "Kraken": os.path.join(GDIR, "Kraken"),
    "CleanSize_vs_BUSCO": os.path.join(GDIR, "CleanSize_vs_BUSCO"),
}
for d in dirs.values():
    os.makedirs(d, exist_ok=True)

cl_csv = os.path.join(dirs["ContigLengths"], f"{SAMPLE}_contigs.csv")
if assembly_path:
    with open(cl_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["contig_id", "length_bp"])
        for rec in SeqIO.parse(assembly_path, "fasta"):
            w.writerow([rec.id, len(rec.seq)])

as_csv = os.path.join(dirs["AssemblyStats"], f"{SAMPLE}_assembly_stats.csv")
with open(as_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    w.writerow(["TotalLength_bp", total_len])
    w.writerow(["N50_bp", n50])
    w.writerow(["MaxContig_bp", max_ctg])
    w.writerow(["NumContigs", n_contigs])

busco_csv = os.path.join(dirs["BUSCO"], f"{SAMPLE}_busco.csv")
with open(busco_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Complete", "Duplicated", "Fragmented", "Missing", "n"])
    w.writerow([busco_c, busco_d, busco_f, busco_m, busco_n])

cov_csv = os.path.join(dirs["Coverage"], f"{SAMPLE}_coverage.csv")
if ddf is not None:
    ddf.to_csv(cov_csv, index=False)

kr_csv = os.path.join(dirs["Kraken"], f"{SAMPLE}_kraken_top10.csv")
if "kdf" in globals() and kdf is not None:
    sdf10 = kdf[kdf["rank"] == "S"].copy().sort_values("perc", ascending=False).head(10)
    sdf10[["name", "perc"]].to_csv(kr_csv, index=False)

csb_csv = os.path.join(dirs["CleanSize_vs_BUSCO"], f"{SAMPLE}_cleanSize_busco.csv")
with open(csb_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Sample", "CleanSize_bp", "BUSCO_Complete"])
    w.writerow([SAMPLE, total_len, busco_c])

lines = []
lines.append("=== Sample summary / TESTOV54 pipeline ===")
lines.append(f"Sample ID          : {SAMPLE}")
lines.append(f"Run time           : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")
lines.append(f"Output directory   : {OUT_DIR}")
lines.append("")
lines.append("== Input & trimming (fastp) ==")
if os.path.exists(fjson):
    lines.append(f"fastp JSON         : {fjson}")
if pairs_before is not None and pairs_after is not None:
    drop = 100.0 * (pairs_before - pairs_after) / max(1.0, pairs_before)
    lines.append(f"Reads (before→after): {int(pairs_before):,} → {int(pairs_after):,}  (-{drop:.1f}%)")
if q30 is not None:
    lines.append(f"Q30 base fraction  : {100.0*float(q30):.1f}%")
if r1len is not None and r2len is not None:
    lines.append(f"Mean read length   : R1={float(r1len):.1f} bp, R2={float(r2len):.1f} bp")

lines.append("")
lines.append("== Assembly (final) ==")
if assembly_path:
    lines.append(f"Assembly used      : {assembly_path}")
if not np.isnan(total_len):
    lines.append(f"Total length       : {total_len/1e6:.2f} Mb")
if not np.isnan(n50):
    lines.append(f"N50                : {n50/1e6:.3f} Mb")
if not np.isnan(max_ctg):
    lines.append(f"Max contig         : {max_ctg/1e6:.3f} Mb")
if not np.isnan(n_contigs):
    lines.append(f"Number of contigs  : {int(n_contigs):,}")

lines.append("")
lines.append("== BUSCO (fungi_odb10) ==")
if not np.isnan(busco_c):
    lines.append(f"Complete           : {busco_c:.1f}%")
    lines.append(f"  └ Duplicated     : {busco_d:.1f}%")
    lines.append(f"Fragmented         : {busco_f:.1f}%")
    lines.append(f"Missing            : {busco_m:.1f}%")
    if not np.isnan(busco_n):
        lines.append(f"Dataset size (n)   : {int(busco_n)}")
else:
    lines.append("BUSCO              : no result found")

lines.append("")
lines.append("== Mapping & coverage (BWA / samtools) ==")
if not np.isnan(depth_mean):
    lines.append(f"Mean coverage      : {depth_mean:.1f}×")
if not np.isnan(depth_median):
    lines.append(f"Median coverage    : {depth_median:.1f}×")
lines.append(f"BAM                : {os.path.join(OUT_DIR, 'mapped.sorted.bam')}")
lines.append(f"Depth table        : {os.path.join(OUT_DIR, 'depth.txt')}")

lines.append("")
lines.append("== Kraken2 top species (rank=S) ==")
if kraken_top:
    for name, perc in kraken_top:
        lines.append(f"{name:30s} {perc:5.1f}%")
else:
    lines.append("No species-level hits in Kraken report (rank=S).")

lines.append("")
lines.append("== Other outputs ==")
lines.append(f"Pilon assembly     : {os.path.join(OUT_DIR, 'pilon_corrected.fasta')}")
lines.append(f"Kraken report      : {os.path.join(OUT_DIR, 'contamination_checks', 'kraken', 'kraken_report.txt')}")
lines.append(f"Sourmash           : {os.path.join(OUT_DIR, 'contamination_checks', 'sourmash', 'assembly.sig')}, {os.path.join(OUT_DIR, 'contamination_checks', 'sourmash', 'gather_results.csv')}")

summary_txt_path = os.path.join(OUT_DIR, "summary.txt")
with open(summary_txt_path, "w", encoding="utf-8") as w:
    w.write("\n".join(lines))

import textwrap
summary_pdf = os.path.join(OUT_DIR, "sample_report.pdf")

wrapped_lines = []
for ln in lines:
    if ln.strip() == "":
        wrapped_lines.append("")
    else:
        wrapped_lines.extend(
            textwrap.wrap(
                ln, width=95,
                break_long_words=False,
                break_on_hyphens=False
            )
        )
text = "\n".join(wrapped_lines)

fig, ax = plt.subplots(figsize=(8.5, 11))
ax.axis("off")
ax.text(
    0.05, 0.97, text,
    va="top", ha="left",
    transform=ax.transAxes,
    fontsize=9.5,
    linespacing=1.35,
    fontstyle="normal",
    fontweight="normal",
    family="monospace"
)
fig.subplots_adjust(left=0.03, right=0.97, top=0.98, bottom=0.03)
fig.savefig(summary_pdf, dpi=300)
plt.close(fig)

print(f"[{SAMPLE}] summary.txt        -> {summary_txt_path}")
print(f"[{SAMPLE}] sample_report.pdf  -> {summary_pdf}")
print(f"[{SAMPLE}] sample_metrics.csv -> {sample_metrics_path}")
print(f"[{SAMPLE}] master_metrics.csv -> {master_path} (append)")
PY_PER_SAMPLE
  deactivate

  touch "$OUT_DIR/.DONE"
done

# =========================
# GLOBAL HEATMAP MATRICES + CONTAMINATION SUMMARY
# =========================
source "$SHARED_VENV/bin/activate"
python3 - << 'PY_HEATMAPS'
import os, glob
import pandas as pd

OUT_ROOT = os.path.expanduser(os.environ["OUT_ROOT"])

def build_matrix(long_files, sample_col, taxa_col, value_col, out_prefix, topN=50):
    if not long_files:
        print(f"[{out_prefix}] No input files found; skipping.")
        return None, None

    long = pd.concat((pd.read_csv(f) for f in long_files), ignore_index=True)
    long[value_col] = pd.to_numeric(long[value_col], errors="coerce").fillna(0.0)

    top_taxa = (
        long.groupby(taxa_col)[value_col]
            .mean()
            .sort_values(ascending=False)
            .head(topN)
            .index
            .tolist()
    )
    long_top = long[long[taxa_col].isin(top_taxa)].copy()

    wide = long_top.pivot_table(
        index=sample_col, columns=taxa_col, values=value_col,
        aggfunc="max", fill_value=0.0
    )

    long_path = os.path.join(OUT_ROOT, f"{out_prefix}_long_all.csv")
    wide_path = os.path.join(OUT_ROOT, f"{out_prefix}_heatmap_matrix_top{topN}.csv")
    long.to_csv(long_path, index=False)
    wide.to_csv(wide_path)

    print(f"[{out_prefix}] wrote: {long_path}")
    print(f"[{out_prefix}] wrote: {wide_path}")
    return long, wide

kraken_files = glob.glob(os.path.join(
    OUT_ROOT, "*_out", "contamination_checks", "kraken", "kraken_species_long.csv"
))
kraken_long, _ = build_matrix(
    kraken_files, "Sample", "name", "perc", "kraken_species", topN=50
)

sm_files = glob.glob(os.path.join(
    OUT_ROOT, "*_out", "contamination_checks", "sourmash", "sourmash_gather_long.csv"
))
sm_long, _ = build_matrix(
    sm_files, "Sample", "match", "value", "sourmash_gather", topN=50
)

rows = []
master_csv = os.path.join(OUT_ROOT, "master_metrics.csv")
master = None
if os.path.exists(master_csv):
    try:
        master = pd.read_csv(master_csv)
    except Exception:
        master = None

samples = sorted({os.path.basename(p).replace("_out","") for p in glob.glob(os.path.join(OUT_ROOT, "*_out"))})

for s in samples:
    row = {"Sample": s}

    k_path = os.path.join(OUT_ROOT, f"{s}_out", "contamination_checks", "kraken", "kraken_species_long.csv")
    if os.path.exists(k_path):
        kdf = pd.read_csv(k_path)
        if not kdf.empty:
            kdf2 = kdf.sort_values("perc", ascending=False)
            row["Top_Kraken_Species"] = str(kdf2.iloc[0]["name"])
            row["Top_Kraken_%Reads"]  = float(kdf2.iloc[0]["perc"])
        else:
            row["Top_Kraken_Species"] = ""
            row["Top_Kraken_%Reads"]  = 0.0
    else:
        row["Top_Kraken_Species"] = ""
        row["Top_Kraken_%Reads"]  = 0.0

    sm_path = os.path.join(OUT_ROOT, f"{s}_out", "contamination_checks", "sourmash", "sourmash_gather_long.csv")
    if os.path.exists(sm_path):
        sdf = pd.read_csv(sm_path)
        if not sdf.empty:
            sdf2 = sdf.sort_values("value", ascending=False)
            row["Top_Sourmash_Match"] = str(sdf2.iloc[0]["match"])
            row["Top_Sourmash_Value"] = float(sdf2.iloc[0]["value"])
        else:
            row["Top_Sourmash_Match"] = ""
            row["Top_Sourmash_Value"] = 0.0
    else:
        row["Top_Sourmash_Match"] = ""
        row["Top_Sourmash_Value"] = 0.0

    if master is not None and "Sample" in master.columns:
        m = master[master["Sample"] == s]
        if not m.empty:
            for col in ["BUSCO_Complete", "TotalLength_bp", "N50_bp", "NumContigs", "Depth_mean", "Fastp_Q30_rate"]:
                if col in m.columns:
                    row[col] = m.iloc[-1][col]

    rows.append(row)

out = pd.DataFrame(rows)
out_path = os.path.join(OUT_ROOT, "contamination_summary.csv")
out.to_csv(out_path, index=False)
print("Wrote contamination summary:", out_path)
PY_HEATMAPS
deactivate || true

fi

# =========================================
# CO-ASSEMBLY (sequential) SECTION
# =========================================

THREADS="${NSLOTS:-8}"

COA_DIR="$OUT_ROOT/newco-assembly"
TOP_DIR="$COA_DIR/top_samples"
SEQ_DIR="$COA_DIR/sequential_assemblies"

mkdir -p "$COA_DIR" "$TOP_DIR" "$SEQ_DIR"

export OUT_ROOT SIF_DIR DB_DIR BUSCO_LINEAGE SM_DB_FUNGI SM_DB_MAMMAL

# =========================
# Step 1: Select Top 12 Samples by BUSCO_Complete
# =========================
source "$SHARED_VENV/bin/activate"
python3 - << 'PY_SELECT'
import os, pandas as pd, shutil, sys

OUT_ROOT = os.environ["OUT_ROOT"]
COA_DIR  = os.path.join(OUT_ROOT, "newco-assembly")
TOP_DIR  = os.path.join(COA_DIR, "top_samples")
CSV      = os.path.join(OUT_ROOT, "master_metrics.csv")

if not os.path.exists(CSV):
    sys.exit(f"ERROR: master_metrics.csv not found: {CSV}")

df = pd.read_csv(CSV)

if "Sample" not in df.columns or "BUSCO_Complete" not in df.columns:
    sys.exit("ERROR: Required columns 'Sample' and 'BUSCO_Complete' missing in master_metrics.csv")

df_nonan = df.dropna(subset=["BUSCO_Complete"])
if df_nonan.empty:
    sys.exit("ERROR: No samples with BUSCO_Complete values found.")

top12 = df_nonan.sort_values("BUSCO_Complete", ascending=False).head(12)

os.makedirs(COA_DIR, exist_ok=True)
os.makedirs(TOP_DIR, exist_ok=True)

rank_txt = os.path.join(COA_DIR, "top12_samples.txt")
rank_csv = os.path.join(COA_DIR, "top12_ranked_samples.csv")

top12["Sample"].to_csv(rank_txt, index=False, header=False)
top12.to_csv(rank_csv, index=False)

print("Top-12 samples by BUSCO_Complete:")
print(top12[["Sample", "BUSCO_Complete"]])

for sample in top12["Sample"]:
    sdir = os.path.join(OUT_ROOT, f"{sample}_out")
    if not os.path.isdir(sdir):
        print(f"WARNING: sample directory not found: {sdir}")
        continue

    dest = os.path.join(TOP_DIR, sample)
    os.makedirs(dest, exist_ok=True)

    candidates_r1 = [
        os.path.join(sdir, "norm_R1.fastq.gz"),
        os.path.join(sdir, "trimmed_R1.fastq.gz"),
        os.path.join(sdir, "clean_R1.fastq.gz"),
        os.path.join(sdir, "raw_outputs", "fastp", "trimmed_R1.fastq.gz"),
    ]
    candidates_r2 = [
        os.path.join(sdir, "norm_R2.fastq.gz"),
        os.path.join(sdir, "trimmed_R2.fastq.gz"),
        os.path.join(sdir, "clean_R2.fastq.gz"),
        os.path.join(sdir, "raw_outputs", "fastp", "trimmed_R2.fastq.gz"),
    ]

    r1 = next((x for x in candidates_r1 if os.path.exists(x)), None)
    r2 = next((x for x in candidates_r2 if os.path.exists(x)), None)

    if not r1 or not r2:
        print(f"WARNING: missing trimmed/clean/normalized reads for {sample}")
        continue

    shutil.copy2(r1, os.path.join(dest, "trimmed_R1.fastq.gz"))
    shutil.copy2(r2, os.path.join(dest, "trimmed_R2.fastq.gz"))

    print(f"Copied reads for {sample}")

PY_SELECT
deactivate

# =========================
# Step 2: Sequential co-assembly function
# =========================
run_seq_coassembly() {
  local K="$1" LABEL="$2"
  local OUT_DIR="$SEQ_DIR/$LABEL"
  mkdir -p "$OUT_DIR"

  local SAMPLES=()
  local i=0
  while IFS= read -r s; do
    [ -z "$s" ] && continue
    SAMPLES+=("$s")
    i=$((i+1))
    [ "$i" -ge "$K" ] && break
  done < "$COA_DIR/top12_samples.txt"

  if [ "${#SAMPLES[@]}" -lt "$K" ]; then
    echo "[coassembly-$LABEL] ERROR: only found ${#SAMPLES[@]} samples in top12 list; need K=$K"
    exit 1
  fi

  local R1_LIST=() R2_LIST=()
  for s in "${SAMPLES[@]}"; do
    local R1="$TOP_DIR/$s/trimmed_R1.fastq.gz"
    local R2="$TOP_DIR/$s/trimmed_R2.fastq.gz"
    [ -f "$R1" ] || { echo "[coassembly-$LABEL] Missing $R1"; exit 1; }
    [ -f "$R2" ] || { echo "[coassembly-$LABEL] Missing $R2"; exit 1; }
    R1_LIST+=("$R1"); R2_LIST+=("$R2")
  done

  zcat "${R1_LIST[@]}" | gzip > "$OUT_DIR/R1_combined.fastq.gz"
  zcat "${R2_LIST[@]}" | gzip > "$OUT_DIR/R2_combined.fastq.gz"

  command -v spades.py >/dev/null 2>&1 || { echo "[coassembly-$LABEL] ERROR: spades.py not found in PATH"; exit 1; }

  spades.py --sc \
    -1 "$OUT_DIR/R1_combined.fastq.gz" \
    -2 "$OUT_DIR/R2_combined.fastq.gz" \
    -o "$OUT_DIR/spades_out" \
    --threads "$THREADS" --memory 64 \
    > "$OUT_DIR/spades.log" 2>&1

  local ASSEMBLY="$OUT_DIR/spades_out/contigs.fasta"
  if [ ! -s "$ASSEMBLY" ]; then
    echo "[coassembly-$LABEL] ERROR: SPAdes did not produce contigs: $ASSEMBLY"
    ls -lh "$OUT_DIR/spades_out" || true
    tail -n 50 "$OUT_DIR/spades.log" || true
    exit 1
  fi

  local KRAKEN_DIR="$OUT_DIR/contamination_checks/kraken"
  mkdir -p "$KRAKEN_DIR"

  singularity exec "$SIF_DIR/kraken2.sif" kraken2 \
    --db "$DB_DIR" \
    --report "$KRAKEN_DIR/kraken_report.txt" \
    --paired "$OUT_DIR/R1_combined.fastq.gz" "$OUT_DIR/R2_combined.fastq.gz" \
    --threads "$THREADS" \
    > "$KRAKEN_DIR/kraken_output.txt" 2>&1

  singularity exec "$SIF_DIR/bwa.sif" bwa index "$ASSEMBLY"
  singularity exec "$SIF_DIR/bwa.sif" bwa mem -t "$THREADS" \
    "$ASSEMBLY" "$OUT_DIR/R1_combined.fastq.gz" "$OUT_DIR/R2_combined.fastq.gz" \
    | singularity exec "$SIF_DIR/samtools.sif" samtools view -bS - \
    > "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools sort -@ "$THREADS" \
    -o "$OUT_DIR/mapped.sorted.bam" "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools index "$OUT_DIR/mapped.sorted.bam"
  singularity exec "$SIF_DIR/samtools.sif" samtools depth -a "$OUT_DIR/mapped.sorted.bam" > "$OUT_DIR/depth.txt"

  singularity exec "$SIF_DIR/pilon.sif" java -Xmx60G -jar /pilon/pilon.jar \
    --genome "$ASSEMBLY" \
    --frags "$OUT_DIR/mapped.sorted.bam" \
    --output "$OUT_DIR/pilon_corrected" \
    --threads "$THREADS" \
    > "$OUT_DIR/pilon.log" 2>&1

  local PILON_ASSEMBLY="$OUT_DIR/pilon_corrected.fasta"
  local FILTERED_ASSEMBLY="$PILON_ASSEMBLY"
  local SEQKIT_SIF="$SIF_DIR/seqkit.sif"

  if [ -f "$SEQKIT_SIF" ]; then
    echo "[coassembly-$LABEL] Filtering contigs < 1000 bp using seqkit..."
    singularity exec "$SEQKIT_SIF" seqkit seq -m 1000 "$PILON_ASSEMBLY" > "$OUT_DIR/pilon_corrected.1kb.fasta"
    FILTERED_ASSEMBLY="$OUT_DIR/pilon_corrected.1kb.fasta"
  else
    echo "[coassembly-$LABEL] Seqkit SIF not found at $SEQKIT_SIF — skipping contig length filter."
  fi

  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_NUM_THREADS=1
  local TS; TS=$(date +%Y%m%d%H%M%S)
  local BUSCO_RUN="busco_${LABEL}_${TS}"
  local BUSCO_LOG="$OUT_DIR/${BUSCO_RUN}.log"

  singularity exec -B "$OUT_DIR:$OUT_DIR" "$SIF_DIR/busco_v6.0.0_cv1.sif" \
    busco -i "$FILTERED_ASSEMBLY" -o "$BUSCO_RUN" -l "$BUSCO_LINEAGE" -m genome --cpu 2 --out_path "$OUT_DIR" \
    > "$BUSCO_LOG" 2>&1

  local SM_OUT="$OUT_DIR/contamination_checks/sourmash"
  mkdir -p "$SM_OUT"
  local SIF_SM="$SIF_DIR/sourmash.sif"

  singularity exec "$SIF_SM" sourmash sketch dna -p k=21,scaled=1000 \
    -o "$SM_OUT/assembly.sig" "$FILTERED_ASSEMBLY"

  [ -f "$SM_DB_FUNGI"  ] || { echo "[coassembly-$LABEL] Missing SM_DB_FUNGI: $SM_DB_FUNGI" >&2; exit 1; }
  [ -f "$SM_DB_MAMMAL" ] || { echo "[coassembly-$LABEL] Missing SM_DB_MAMMAL: $SM_DB_MAMMAL" >&2; exit 1; }

  singularity exec "$SIF_SM" sourmash gather \
    "$SM_OUT/assembly.sig" \
    "$SM_DB_FUNGI" \
    "$SM_DB_MAMMAL" \
    -o "$SM_OUT/gather_results.csv"

  echo "[coassembly-$LABEL] Finished co-assembly K=$K (OUT_DIR=$OUT_DIR)"
}

# =========================
# Step 3: Run sequential co-assemblies (2→12)
# =========================
run_seq_coassembly 2 double
run_seq_coassembly 3 triple
run_seq_coassembly 4 quadruple
run_seq_coassembly 5 pentuple
run_seq_coassembly 6 sextuple
run_seq_coassembly 7 septuple
run_seq_coassembly 8 octuple
run_seq_coassembly 9 nonuple
run_seq_coassembly 10 decuple
run_seq_coassembly 11 elevenfold
run_seq_coassembly 12 twelvefold

# =========================
# Step 4: BUSCO completeness trend plot
# =========================
source "$SHARED_VENV/bin/activate"
python3 - << 'PY_TREND'
import os, re, glob
import numpy as np
import matplotlib.pyplot as plt

OUT_ROOT = os.path.expanduser(os.environ.get("OUT_ROOT", "results/andeep_250k"))
COA_DIR  = os.path.join(OUT_ROOT, "newco-assembly")
SEQ_DIR  = os.path.join(COA_DIR, "sequential_assemblies")

labels = [("double",2),("triple",3),("quadruple",4),("pentuple",5),("sextuple",6),("septuple",7),("octuple",8),("nonuple",9),("decuple",10),("elevenfold",11),("twelvefold",12)]

x = []
complete = []
fragmented = []

for lbl,k in labels:
    outdir = os.path.join(SEQ_DIR,lbl)
    logs = sorted(
        glob.glob(os.path.join(outdir,"logs","busco_*.log")) +
        glob.glob(os.path.join(outdir,"busco_*.log")),
        key=os.path.getmtime, reverse=True
    )
    comp, frag = np.nan, np.nan
    for lg in logs:
        for line in open(lg, encoding='utf-8', errors='ignore'):
            m = re.search(r"C:(\d+\.\d+)%.*F:(\d+\.\d+)%", line)
            if m:
                comp, frag = float(m.group(1)), float(m.group(2))
                break
        if not np.isnan(comp):
            break
    x.append(k)
    complete.append(comp)
    fragmented.append(frag)

blue = "#0072BC"
yellow = "#FFC740"

plt.figure(figsize=(8,5))
plt.scatter(x, complete, color=blue, label='Complete', s=100, marker='o')
plt.plot(x, complete, color=blue, linestyle='--', alpha=0.7)
for xi, yi in zip(x, complete):
    if not np.isnan(yi):
        plt.text(xi, yi + 1, f"{yi:.1f}%", ha='center', va='bottom', fontsize=9, color='black')

plt.scatter(x, fragmented, color=yellow, label='Fragmented', s=100, marker='^')
plt.plot(x, fragmented, color=yellow, linestyle='--', alpha=0.7)
for xi, yi in zip(x, fragmented):
    if not np.isnan(yi):
        plt.text(xi, yi - 2, f"{yi:.1f}%", ha='center', va='top', fontsize=9, color='black')

plt.xticks(x, ["Double","Triple","Quadruple","Pentuple","Sextuple","Septuple","Octuple","Nonuple","Decuple","Elevenfold","Twelvefold"])
plt.xlabel("Number of SAGs in co-assembly")
plt.ylabel("BUSCO %")
plt.title("BUSCO completeness across co-assemblies")
plt.ylim(0, 105)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()

out_png = os.path.join(COA_DIR, "sequential_results.png")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved BUSCO line plot: {out_png}")
PY_TREND
deactivate || true

echo "Sequential co-assembly workflow complete."
