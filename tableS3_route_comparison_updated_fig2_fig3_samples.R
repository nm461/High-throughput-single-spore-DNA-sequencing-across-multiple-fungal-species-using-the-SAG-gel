#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

data_root <- Sys.getenv(
  "PAPER_DATA_ROOT",
  unset = normalizePath(file.path(getwd(), "PaperDataFiles"), mustWork = FALSE)
)


repo_root <- file.path(data_root, "graphs", "standardisedpaperfigures")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) repo_root <- normalizePath(args[1], mustWork = TRUE)

outdir <- file.path(repo_root, "outputs", "tableS3_route_comparison_updated_fig2_fig3_samples")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_sample_col <- function(path) {
  if (!file.exists(path)) stop("Missing sample-set file: ", path)

  read_csv(path, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    transmute(sample = as.character(sample)) %>%
    distinct()
}

read_nn_metrics <- function(path, species, route) {
  if (!file.exists(path)) stop("Missing assembly metrics file: ", path)

  ext <- tools::file_ext(path)
  df <- if (ext == "xlsx") readxl::read_xlsx(path) else read_csv(path, show_col_types = FALSE)

  df %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    mutate(
      sample = as.character(sample),
      species = species,
      route = route,
      busco_complete = as.numeric(busco_complete),
      n50_kbp = as.numeric(n50_bp) / 1000,
      contig_count = as.numeric(numcontigs),
      reads_after = as.numeric(reads_after),
      self_mapping_mean_depth = as.numeric(depth_mean),
      total_length_mbp = as.numeric(totallength_bp) / 1e6
    ) %>%
    select(
      sample, species, route,
      busco_complete, n50_kbp, contig_count, reads_after,
      self_mapping_mean_depth, total_length_mbp
    )
}

inputs <- tibble::tribble(
  ~species, ~route, ~assembly_metrics, ~refmap_metrics, ~fig4_sample_set,
  "A. niger", "direct",
  file.path(data_root, "assembly", "single_SAG_nonnorm", "AN_master_metrics.xlsx"),
  file.path(data_root, "reference_mapping", "refmap_metrics", "AN_shallow_direct_refmap_metrics_200k.csv"),
  file.path(data_root, "reference_mapping", "bin_qc", "AN_shallow_direct_bin_qc_10kb.csv"),
  "A. niger", "indirect",
  file.path(data_root, "assembly", "single_SAG_nonnorm", "ANdeep_master_metrics.csv"),
  file.path(data_root, "reference_mapping", "refmap_metrics", "AN_deep_indirect_refmap_metrics_200k.csv"),
  file.path(data_root, "reference_mapping", "bin_qc", "AN_deep_indirect_bin_qc_10kb.csv"),
  "C. nagasakiense", "direct",
  file.path(data_root, "assembly", "single_SAG_nonnorm", "CNTW_master_metrics.xlsx"),
  file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_shallow_direct_refmap_metrics_200k.csv"),
  file.path(data_root, "reference_mapping", "bin_qc", "CNTW_shallow_direct_bin_qc_10kb.csv"),
  "C. nagasakiense", "indirect",
  file.path(data_root, "assembly", "single_SAG_nonnorm", "CNTWdeep_master_metrics.csv"),
  file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_deep_indirect_refmap_metrics_200k.csv"),
  file.path(data_root, "reference_mapping", "bin_qc", "CNTW_deep_indirect_bin_qc_10kb.csv")
)

filtered_inputs <- bind_rows(lapply(seq_len(nrow(inputs)), function(i) {
  row <- inputs[i, ]

  matched_samples <- read_sample_col(row$fig4_sample_set) %>%
    inner_join(read_sample_col(row$refmap_metrics), by = "sample") %>%
    filter(sample != "AN_7")

  read_nn_metrics(row$assembly_metrics, row$species, row$route) %>%
    filter(sample != "AN_7") %>%
    semi_join(matched_samples, by = "sample")
}))

expected_counts <- tibble::tribble(
  ~species, ~route, ~expected_n,
  "A. niger", "direct", 175L,
  "A. niger", "indirect", 10L,
  "C. nagasakiense", "direct", 170L,
  "C. nagasakiense", "indirect", 15L
)

count_check <- filtered_inputs %>%
  distinct(sample, species, route) %>%
  count(species, route, name = "n") %>%
  right_join(expected_counts, by = c("species", "route")) %>%
  mutate(n = coalesce(n, 0L))

if (any(count_check$n != count_check$expected_n)) {
  stop(
    "Sample count mismatch:\n",
    paste(
      sprintf(
        "%s / %s expected %d found %d",
        count_check$species, count_check$route, count_check$expected_n, count_check$n
      ),
      collapse = "\n"
    )
  )
}

metrics <- tibble::tribble(
  ~Metric, ~column, ~unit,
  "BUSCO completeness", "busco_complete", "%",
  "N50", "n50_kbp", "kbp",
  "Number of contigs", "contig_count", "",
  "Reads after fastp", "reads_after", "",
  "Self-mapping mean depth", "self_mapping_mean_depth", "x",
  "Total assembled length", "total_length_mbp", "Mbp"
)

format_median <- function(x, unit) {
  if (unit == "%") {
    sprintf("%.1f%%", x)
  } else if (unit == "kbp") {
    sprintf("%.2f kbp", x)
  } else if (unit == "Mbp") {
    sprintf("%.2f Mbp", x)
  } else if (unit == "x") {
    sprintf("%.2fx", x)
  } else {
    format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
}

table_s3 <- bind_rows(lapply(c("A. niger", "C. nagasakiense"), function(sp) {
  bind_rows(lapply(seq_len(nrow(metrics)), function(i) {
    metric <- metrics[i, ]
    species_dat <- filtered_inputs %>% filter(species == sp)

    direct_values <- species_dat %>%
      filter(route == "direct") %>%
      pull(.data[[metric$column]])

    indirect_values <- species_dat %>%
      filter(route == "indirect") %>%
      pull(.data[[metric$column]])

    direct_values <- direct_values[is.finite(direct_values)]
    indirect_values <- indirect_values[is.finite(indirect_values)]

    test <- wilcox.test(direct_values, indirect_values, alternative = "two.sided", exact = FALSE)

    tibble(
      Species = sp,
      Metric = metric$Metric,
      `n (direct)` = length(direct_values),
      `n (indirect)` = length(indirect_values),
      `Median (direct)` = format_median(median(direct_values), metric$unit),
      `Median (indirect)` = format_median(median(indirect_values), metric$unit),
      U = as.numeric(test$statistic),
      `p (raw)` = test$p.value
    )
  }))
})) %>%
  mutate(
    `p (BH-FDR)` = p.adjust(`p (raw)`, method = "BH"),
    Significant = ifelse(`p (BH-FDR)` < 0.05, "Yes", "No")
  )

write_csv(table_s3, file.path(outdir, "Table_S3_route_comparison_updated.csv"))
write_csv(filtered_inputs, file.path(outdir, "Table_S3_filtered_input_samples.csv"))

print(table_s3, n = Inf)
message("Wrote: ", normalizePath(file.path(outdir, "Table_S3_route_comparison_updated.csv")))
message("Wrote: ", normalizePath(file.path(outdir, "Table_S3_filtered_input_samples.csv")))
