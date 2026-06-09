# Core Pipeline Scripts

These scripts document the core computational workflows used for the SAG-gel single-spore fungal genomics study. They are intended as reproducibility scripts for manuscript review/publication, not as polished general-purpose software.

The scripts were originally run on an SGE cluster with Apptainer/Singularity containers. Personal absolute input/output paths have been replaced with environment variables so the same commands can be rerun on another system using deposited FASTQ files, reference assemblies, and database resources.

## Scripts

| Script | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| `Spore_assembly.sh` | Direct-route shallow single-SAG assembly workflow. Runs fastp, optional host decontamination, SPAdes single-cell assembly, BWA/SAMtools self-mapping, Pilon polishing, contig filtering, BUSCO, and per-sample assembly metrics. | Paired FASTQ files named `spore_1_sc-*_1.fastq.gz` and `spore_1_sc-*_2.fastq.gz` | Per-sample assembly folders, polished assemblies, BUSCO logs, `sample_metrics.csv`, `master_metrics.csv` |
| `Normalisedassembly.sh` | Indirect-route/read-normalized single-SAG assembly workflow. Runs fastp, optional host decontamination, 250k read-pair downsampling with rasusa, SPAdes, contamination screening, polishing, BUSCO, and summary metrics. | Paired FASTQ files under per-sample directories, e.g. `AN*/sample_R1.fastq.gz` and `AN*/sample_R2.fastq.gz` | Per-sample assemblies, BUSCO and contamination summaries, `master_metrics.csv`, co-assembly summaries when enabled |
| `Referencemapping.sh` | Reference mapping and coverage-uniformity workflow. Runs BWA mapping, SAMtools metrics, read-count normalization, Picard metrics, mosdepth 10 kb/50 kb bin depth, GC/bin QC, Gini metrics, and Lorenz-curve tables. | Per-sample paired FASTQ directories, reference FASTA, comparator `master_refmap_metrics.csv` for normalization target | `master_refmap_metrics.csv`, `normalization_plan.csv`, `master_bin_qc_metrics_*`, `master_lorenz_curve_*`, `master_picard_metrics.csv` |
| `Lorenzdatageneration.sh` | Utility script to regenerate Lorenz-curve CSV files from existing mosdepth region outputs. | Sample directories containing `mosdepth_10kb.regions.bed.gz` and/or `mosdepth_50kb.regions.bed.gz` | Per-sample and master Lorenz-curve CSVs |
| `Bootstraprandom_coasssembly.sh` | Random bootstrap co-assembly workflow for direct-route shallow SAGs. Selects random eligible SAG subsets, concatenates reads, runs SPAdes co-assembly, Pilon, contig filtering, BUSCO, and records BUSCO completeness per subset size. | Completed single-SAG assembly output tree containing `master_metrics.csv` and per-sample normalized/clean FASTQs | Iteration-level subset manifests, per-subset BUSCO logs, `busco_metrics.csv` |

Software versions and container-version check commands are listed in `software_versions.txt`.

## Quick Start

Create or edit an environment file for the compute system:

```bash
cp env.example env.local
vim env.local
source env.local
```

Run a script directly:

```bash
bash Spore_assembly.sh
```

or submit it to an SGE cluster:

```bash
mkdir -p logs
qsub Spore_assembly.sh
```

For array jobs, set the array range in the SGE header or at submission time. For example:

```bash
mkdir -p logs
qsub -t 1-193 Spore_assembly.sh
qsub -t 1-100 Bootstraprandom_coasssembly.sh
```

## Running Multiple Datasets

The scripts are designed to be rerun with different `DATA_DIR`, `OUT_ROOT`, `MAP_ROOT`, `SAMPLES_ROOT`, and `REF_SOURCE` values. For example:

```bash
source env.local

DATA_DIR="$PWD/data/direct_aniger_fastq" \
OUT_ROOT="$PWD/results/direct_aniger_assembly" \
MAX_SAMPLES=193 \
bash Spore_assembly.sh

DATA_DIR="$PWD/data/indirect_aniger_fastq" \
OUT_ROOT="$PWD/results/indirect_aniger_250k" \
NORMALIZE_READS=250000 \
MIN_READS=250000 \
MAX_SAMPLES=16 \
bash Normalisedassembly.sh

SAMPLES_ROOT="$PWD/data/direct_aniger_fastq_by_sample" \
MAP_ROOT="$PWD/results/direct_aniger_refmap" \
REF_SOURCE="$PWD/references/GCF_000002855.4_ASM285v2_genomic.fna" \
REFERENCE_MASTER="$PWD/results/direct_cntw_refmap/master_refmap_metrics.csv" \
bash Referencemapping.sh
```

This pattern is preferable to editing paths inside the scripts for each dataset.

## Expected Input Layout

`Spore_assembly.sh` expects paired FASTQ files directly in `DATA_DIR`:

```text
DATA_DIR/
  spore_1_sc-001_1.fastq.gz
  spore_1_sc-001_2.fastq.gz
  spore_1_sc-002_1.fastq.gz
  spore_1_sc-002_2.fastq.gz
```

`Normalisedassembly.sh` expects per-sample subdirectories:

```text
DATA_DIR/
  AN_001/
    AN_001_R1.fastq.gz
    AN_001_R2.fastq.gz
  AN_002/
    AN_002_R1.fastq.gz
    AN_002_R2.fastq.gz
```

`Referencemapping.sh` expects per-sample subdirectories where the directory name is also the sample prefix:

```text
SAMPLES_ROOT/
  sample_sc-001/
    sample_sc-001_R1.fastq.gz
    sample_sc-001_R2.fastq.gz
```

`Lorenzdatageneration.sh` expects existing reference-mapping sample output directories:

```text
MAP_ROOT/
  sample_sc-001/
    mosdepth_10kb.regions.bed.gz
    mosdepth_50kb.regions.bed.gz
```

`Bootstraprandom_coasssembly.sh` expects the output structure produced by the single-SAG assembly workflow:

```text
OUT_ROOT/
  master_metrics.csv
  sample_001_out/
    .DONE
    norm_R1.fastq.gz
    norm_R2.fastq.gz
```

## Key Parameters

| Variable | Meaning |
|---|---|
| `SIF_DIR` | Directory containing Apptainer/Singularity images such as `fastp.sif`, `bwa.sif`, `samtools.sif`, `pilon.sif`, `busco_v6.0.0_cv1.sif`, `seqkit.sif`, `rasusa.sif`, `kraken2.sif`, and `sourmash.sif`. |
| `SPADES_BIN` | Directory containing `spades.py` from SPAdes v4.2.0. |
| `BUSCO_LINEAGE` | BUSCO lineage directory, here `fungi_odb10`. |
| `HOST_REF` | Optional host/human masked FASTA for read-level decontamination. If absent, host decontamination is skipped. |
| `NORMALIZE_READS` | Target read-pair count for rasusa downsampling in `Normalisedassembly.sh`. Default: `250000`. |
| `REF_SOURCE` | Reference FASTA used by `Referencemapping.sh`. |
| `REFERENCE_MASTER` | Comparator mapping table used to define the normalization target for reference-mapping comparisons. |
| `SUBSET_SIZE`, `MIN_SAG`, `MAX_SAG`, `BASE_SEED` | Co-assembly bootstrap subset settings. |

## Notes for Reproducibility

The scripts intentionally keep the exact workflow order, tool choices, and major thresholds used in the manuscript analyses. They use environment variables for paths so that deposited datasets can be substituted without modifying script internals.

Some resources are large and are not included in this script bundle, including reference FASTA files, raw FASTQ files, Kraken2 databases, Sourmash databases, BUSCO lineage files, and Apptainer/Singularity images. These should be cited or described separately in the manuscript data availability statement.

The plotting scripts and final figure rendering are separate from this core-pipeline bundle. This folder covers the upstream processing described as read processing, assembly, reference mapping, coverage profiling, and co-assembly.
