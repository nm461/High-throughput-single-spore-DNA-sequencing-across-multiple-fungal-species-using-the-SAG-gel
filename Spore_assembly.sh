#!/bin/bash
#$ -N SPOREAssembly
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=400:00:00
#$ -l s_vmem=16G,mem_req=16G
#$ -pe def_slot 8
#$ -o logs/SPOREAssembly.out
#$ -e logs/SPOREAssembly.log
#$ -t 1-193           # one task per sample

set -x

START_AT="${START_AT:-samples}"   # all | coassembly | samples

# =========================
# Modules
# =========================
module use /usr/local/package/modulefiles
module load apptainer
module load blast+/2.15.0
module load repeatmasker/4.1.6 || true

export PATH="$HOME/tools/rmblast/rmblast-2.14.1/bin:$PATH"

# =========================
# Tools on PATH (SPAdes)
# =========================
SPADES_BIN="${SPADES_BIN:-$HOME/bin/SPAdes-4.2.0-Linux/bin}"
export PATH="$SPADES_BIN:$PATH"

command -v spades.py >/dev/null 2>&1 || { echo "ERROR: spades.py not found in PATH ($SPADES_BIN)"; exit 1; }

# =========================
# Directories & inputs
# =========================
OUT_ROOT="${OUT_ROOT:-$PWD/results/spore_assembly}"
SIF_DIR="${SIF_DIR:-$HOME/sif_images}"
DATA_DIR="${DATA_DIR:-$PWD/data/spore_assembly_fastq}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-$HOME/fungi_odb10}"
HOST_REF="${HOST_REF:-$HOME/databases/human_masked/human_masked.fasta}"

MAX_SAMPLES="${MAX_SAMPLES:-193}"

mkdir -p "$OUT_ROOT/logs"
mkdir -p "$OUT_ROOT"
export SIF_DIR
export OUT_ROOT

# =========================
# Sample list (for array-job indexing)
# =========================
SAMPLE_LIST="$OUT_ROOT/sample_list.txt"
if [ ! -f "$SAMPLE_LIST" ]; then
  (
    flock -x 200
    if [ ! -f "$SAMPLE_LIST" ]; then
      shopt -s nullglob
      i=0
      for _R1 in "$DATA_DIR"/spore_1_sc-*_1.fastq.gz; do
        [ "$i" -ge "$MAX_SAMPLES" ] && break
        echo "$_R1"
        i=$((i+1))
      done > "$SAMPLE_LIST"
      echo "Generated sample_list.txt: $(wc -l < "$SAMPLE_LIST") samples → $SAMPLE_LIST"
    fi
  ) 200>"$OUT_ROOT/.sample_list.lock"
fi

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

if [ -n "${SGE_TASK_ID:-}" ]; then
  _array_r1=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
  if [ -z "$_array_r1" ] || [ ! -f "$_array_r1" ]; then
    echo "Task $SGE_TASK_ID: no file at line $SGE_TASK_ID of sample_list.txt — nothing to do."
    exit 0
  fi
  _R1_LIST=("$_array_r1")
else
  _R1_LIST=("$DATA_DIR"/spore_1_sc-*_1.fastq.gz)
fi

for R1 in "${_R1_LIST[@]}"; do
  if [ "$processed" -ge "$MAX_SAMPLES" ]; then
    echo "Reached MAX_SAMPLES=$MAX_SAMPLES; stopping."
    break
  fi

  sample=$(basename "$R1" _1.fastq.gz)
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"

  if [ ! -f "$R2" ]; then
    echo "WARNING: missing R2 for sample $sample — skipping."
    continue
  fi

  processed=$((processed + 1))

  OUT_DIR="$OUT_ROOT/${sample}_out"
  mkdir -p "$OUT_DIR"

  if [ -f "$OUT_DIR/.DONE" ]; then
    echo "==== $sample already done; skipping ===="
    continue
  fi
  if [ -f "$OUT_DIR/.FAILED" ]; then
    echo "==== $sample previously failed; skipping (remove $OUT_DIR to retry) ===="
    continue
  fi
  echo "==== Processing $sample ===="
  export sample R1 R2 OUT_DIR
  (
  set -eo pipefail

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

  R1_USE="$R1_CLEAN"
  R2_USE="$R2_CLEAN"

  # ---------- 2. SPAdes (single-cell mode) ----------
  export PATH="$SPADES_BIN:$PATH"
  ASSEMBLY="$OUT_DIR/spades_out/contigs.fasta"

  if [ -f "$ASSEMBLY" ] && grep -qm1 "^>" "$ASSEMBLY"; then
    echo "[$sample] SPAdes contigs already exist and are valid — skipping SPAdes"
  else
    rm -rf "$OUT_DIR/spades_out"
    spades.py \
      --sc \
      -1 "$R1_USE" \
      -2 "$R2_USE" \
      -o "$OUT_DIR/spades_out" \
      --threads "$NSLOTS" \
      --memory 120 \
      > "$OUT_DIR/spades.log" 2>&1
  fi

  if [ ! -s "$ASSEMBLY" ] || ! grep -qm1 "^>" "$ASSEMBLY"; then
    echo "ERROR [$sample]: SPAdes produced no contigs — skipping sample"
    exit 1
  fi

  # ---------- 3. BWA/Samtools ----------
  singularity exec "$SIF_DIR/bwa.sif" bwa index "$ASSEMBLY"
  singularity exec "$SIF_DIR/bwa.sif" bwa mem -t "$NSLOTS" \
    "$ASSEMBLY" "$R1_USE" "$R2_USE" \
    | singularity exec "$SIF_DIR/samtools.sif" samtools view -bS - \
    > "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools sort -@ "$NSLOTS" \
    -o "$OUT_DIR/mapped.sorted.bam" "$OUT_DIR/mapped.bam"

  singularity exec "$SIF_DIR/samtools.sif" samtools index "$OUT_DIR/mapped.sorted.bam"

  # ---------- 4. Pilon ----------
  singularity exec "$SIF_DIR/pilon.sif" java -Xmx120G -jar /pilon/pilon.jar \
    --genome "$ASSEMBLY" \
    --frags "$OUT_DIR/mapped.sorted.bam" \
    --output "$OUT_DIR/pilon_corrected" \
    --threads "$NSLOTS" \
    > "$OUT_DIR/pilon.log" 2>&1

  PILON_ASSEMBLY="$OUT_DIR/pilon_corrected.fasta"
  if [ ! -s "$PILON_ASSEMBLY" ]; then
    echo "ERROR [$sample]: pilon produced no output — skipping sample"
    exit 1
  fi

  # ---------- 4.5 Contig filtering (>= 1000 bp) ----------
  FILTERED_ASSEMBLY="$PILON_ASSEMBLY"
  SEQKIT_SIF="$SIF_DIR/seqkit.sif"

  if [ -f "$SEQKIT_SIF" ]; then
    echo "Filtering contigs < 1000 bp using seqkit (via SIF)..."
    singularity exec "$SEQKIT_SIF" seqkit seq -m 1000 "$PILON_ASSEMBLY" > "$OUT_DIR/pilon_corrected.1kb.fasta"
    FILTERED_ASSEMBLY="$OUT_DIR/pilon_corrected.1kb.fasta"
  else
    echo "Seqkit SIF not found at $SEQKIT_SIF — skipping contig length filter."
  fi

  SEQKIT_FALLBACK=no
  if [ ! -s "$FILTERED_ASSEMBLY" ] || ! grep -qm1 "^[ACGTNacgtn]" "$FILTERED_ASSEMBLY"; then
    echo "WARNING [$sample]: seqkit filtered all contigs; falling back to pilon_corrected.fasta"
    FILTERED_ASSEMBLY="$PILON_ASSEMBLY"
    SEQKIT_FALLBACK=yes
  fi

  if ! grep -qm1 "^[ACGTNacgtn]" "$FILTERED_ASSEMBLY"; then
    echo "ERROR [$sample]: assembly has no nucleotide sequences; cannot run BUSCO — skipping sample"
    exit 1
  fi

  # ---------- 5. BUSCO ----------
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
    echo "ASSEMBLY: $FILTERED_ASSEMBLY"; ls -lh "$FILTERED_ASSEMBLY" || true
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

  # ---------- 6. Per-sample report + CSV metrics ----------
  source "$SHARED_VENV/bin/activate"
  export OUT_DIR OUT_ROOT SIF_DIR sample SEQKIT_FALLBACK
  python3 - << 'PY_PER_SAMPLE'
import os, json, re, glob, datetime, csv, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
ddf = None

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
    "Seqkit_fallback",
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
    "Seqkit_fallback": os.environ.get("SEQKIT_FALLBACK", "no"),
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
    "ContigLengths": os.path.join(GDIR, "ContigLengths"),
    "AssemblyStats": os.path.join(GDIR, "AssemblyStats"),
    "BUSCO": os.path.join(GDIR, "BUSCO"),
    "Coverage": os.path.join(GDIR, "Coverage"),
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

import textwrap
lines = []
lines.append("=== Sample summary / spore1 pipeline ===")
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
lines.append("== Other outputs ==")
lines.append(f"Pilon assembly     : {os.path.join(OUT_DIR, 'pilon_corrected.fasta')}")

summary_txt_path = os.path.join(OUT_DIR, "summary.txt")
with open(summary_txt_path, "w", encoding="utf-8") as w:
    w.write("\n".join(lines))

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
  ) # end per-sample subshell
  SAMPLE_EXIT=$?
  if [ "$SAMPLE_EXIT" -ne 0 ]; then
    echo "WARNING: [$sample] failed (exit $SAMPLE_EXIT) — recording N/A and marking .FAILED"
    touch "$OUT_DIR/.FAILED"
    source "$SHARED_VENV/bin/activate"
    python3 - << 'PY_NA_RECORD'
import csv, os
OUT_ROOT = os.environ.get("OUT_ROOT", "")
sample   = os.environ.get("sample", "")
path     = os.path.join(OUT_ROOT, "master_metrics.csv")
cols = ["Sample","Assembly_path","TotalLength_bp","N50_bp","MaxContig_bp","NumContigs",
        "BUSCO_Complete","BUSCO_Duplicated","BUSCO_Fragmented","BUSCO_Missing","BUSCO_n",
        "Depth_mean","Depth_median","Fastp_Q30_rate","Fastp_R1_mean_len","Fastp_R2_mean_len",
        "Reads_before","Reads_after","Seqkit_fallback"]
row = {c: "N/A" for c in cols}
row["Sample"] = sample
write_header = not os.path.exists(path)
with open(path, "a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    if write_header:
        w.writeheader()
    w.writerow(row)
print(f"N/A record written for {sample}")
PY_NA_RECORD
    deactivate || true
  fi
done

fi

if [ "$START_AT" = "samples" ]; then
  echo "START_AT=samples — skipping co-assembly."
  exit 0
fi
