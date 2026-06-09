#!/bin/bash
#$ -N CNTWdeep_lorenz_csvs
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l s_vmem=4G,mem_req=4G
#$ -pe def_slot 1
#$ -o logs/CNTW_lorenz_csvs.out
#$ -e logs/CNTW_lorenz_csvs.err

set -e
set -x
set -o pipefail

# ============================================================
# CONFIG
# ============================================================

MAP_ROOT="${MAP_ROOT:-$PWD/results/reference_mapping}"
LOG_ROOT="${LOG_ROOT:-$MAP_ROOT/logs}"

mkdir -p "$MAP_ROOT" "$LOG_ROOT"

# ============================================================
# STEP: build Lorenz-curve CSVs from mosdepth region outputs
# ============================================================

export MAP_ROOT

python3 - <<'PY'
import os
import csv
import gzip
from pathlib import Path

MAP_ROOT = Path(os.environ["MAP_ROOT"])

BIN_LABELS = ["10kb", "50kb"]
INPUT_SUFFIX = {
    "10kb": "mosdepth_10kb.regions.bed.gz",
    "50kb": "mosdepth_50kb.regions.bed.gz",
}

sample_dirs = sorted([p for p in MAP_ROOT.iterdir() if p.is_dir() and p.name != "reference"])

if not sample_dirs:
    raise SystemExit(f"ERROR: no sample directories found under {MAP_ROOT}")

for bin_label in BIN_LABELS:
    master_out = MAP_ROOT / f"master_lorenz_curve_{bin_label}.csv"

    header = [
        "Sample",
        "BinLabel",
        "BinRank",
        "BinsTotal",
        "BinFraction",
        "CumulativeBinFraction",
        "Depth",
        "CumulativeDepth",
        "CumulativeDepthFraction",
    ]

    master_rows = []

    for sd in sample_dirs:
        sample = sd.name
        in_file = sd / INPUT_SUFFIX[bin_label]

        if not in_file.exists():
            print(f"WARNING: missing {in_file}; skipping {sample} for {bin_label}")
            continue

        depths = []

        with gzip.open(in_file, "rt") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue

                try:
                    depth = float(parts[3])
                except ValueError:
                    continue

                depths.append(depth)

        if not depths:
            print(f"WARNING: no usable depths found in {in_file}; skipping {sample} for {bin_label}")
            continue

        depths.sort()  # ascending for Lorenz curve
        n = len(depths)
        total_depth = sum(depths)

        out_file = sd / f"lorenz_curve_{bin_label}.csv"

        rows = []
        cumulative_depth = 0.0

        for i, depth in enumerate(depths, start=1):
            cumulative_depth += depth
            bin_fraction = 1.0 / n
            cumulative_bin_fraction = i / n
            cumulative_depth_fraction = (cumulative_depth / total_depth) if total_depth > 0 else 0.0

            row = [
                sample,
                bin_label,
                i,
                n,
                f"{bin_fraction:.10f}",
                f"{cumulative_bin_fraction:.10f}",
                f"{depth:.10f}",
                f"{cumulative_depth:.10f}",
                f"{cumulative_depth_fraction:.10f}",
            ]
            rows.append(row)
            master_rows.append(row)

        with out_file.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)

        print(f"Wrote per-sample Lorenz CSV: {out_file}")

    with master_out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(master_rows)

    print(f"Wrote master Lorenz CSV: {master_out}")

print("DONE: generated Lorenz curve CSVs from mosdepth region outputs")
PY

echo "Lorenz CSV generation complete."
echo "10kb master: $MAP_ROOT/master_lorenz_curve_10kb.csv"
echo "50kb master: $MAP_ROOT/master_lorenz_curve_50kb.csv"
