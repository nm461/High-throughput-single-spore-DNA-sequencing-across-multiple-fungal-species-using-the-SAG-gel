#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

data_root <- Sys.getenv(
  "PAPER_DATA_ROOT",
  unset = normalizePath(file.path(getwd(), "PaperDataFiles"), mustWork = FALSE)
)


input_table <- file.path(data_root, "supplementary_tables", "master_sample_table.tsv")
repo_root <- file.path(data_root, "graphs", "standardisedpaperfigures")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) input_table <- normalizePath(args[1], mustWork = TRUE)
if (length(args) >= 2) repo_root <- normalizePath(args[2], mustWork = TRUE)

outdir <- file.path(repo_root, "outputs", "table_sample_qc_split")
if (length(args) >= 3) outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_sample_id <- function(x) {
  str_trim(str_remove_all(as.character(x), "\\[.*?\\]"))
}

read_sample_set <- function(path) {
  if (!file.exists(path)) stop("Missing sample-set file: ", path)

  read_csv(path, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    transmute(sample = clean_sample_id(sample)) %>%
    distinct()
}

sample_set_files <- tibble::tribble(
  ~Species, ~Route, ~sample_set_file, ~expected_n,
  "Aspergillus niger", "Direct",
  file.path(data_root, "reference_mapping", "bin_qc", "AN_shallow_direct_bin_qc_10kb.csv"),
  175L,
  "Aspergillus niger", "Indirect",
  file.path(data_root, "reference_mapping", "bin_qc", "AN_deep_indirect_bin_qc_10kb.csv"),
  10L,
  "Colletotrichum nagasakiense", "Direct",
  file.path(data_root, "reference_mapping", "bin_qc", "CNTW_shallow_direct_bin_qc_10kb.csv"),
  170L,
  "Colletotrichum nagasakiense", "Indirect",
  file.path(data_root, "reference_mapping", "bin_qc", "CNTW_deep_indirect_bin_qc_10kb.csv"),
  15L
)

qc_sample_set <- bind_rows(lapply(seq_len(nrow(sample_set_files)), function(i) {
  row <- sample_set_files[i, ]

  read_sample_set(row$sample_set_file) %>%
    filter(sample != "AN_7") %>%
    mutate(Species = row$Species, Route = row$Route)
}))

count_check <- qc_sample_set %>%
  count(Species, Route, name = "n") %>%
  right_join(sample_set_files %>% select(Species, Route, expected_n), by = c("Species", "Route")) %>%
  mutate(n = coalesce(n, 0L))

if (any(count_check$n != count_check$expected_n)) {
  stop(
    "QC sample-set count mismatch:\n",
    paste(
      sprintf(
        "%s / %s expected %d found %d",
        count_check$Species, count_check$Route, count_check$expected_n, count_check$n
      ),
      collapse = "\n"
    )
  )
}

master_table <- read_tsv(input_table, show_col_types = FALSE) %>%
  rename_with(~ str_replace_all(.x, fixed("×"), "x")) %>%
  rename_with(~ str_replace_all(.x, fixed("Ã—"), "x")) %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, fixed("×"), "x"))) %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, fixed("Ã—"), "x"))) %>%
  mutate(`Matched sample ID` = clean_sample_id(`Sample ID`)) %>%
  relocate(`Matched sample ID`, .after = `Sample ID`)

all_with_status <- master_table %>%
  left_join(
    qc_sample_set %>%
      transmute(Species, Route, `Matched sample ID` = sample, qc_passed = TRUE),
    by = c("Species", "Route", "Matched sample ID")
  ) %>%
  mutate(
    qc_passed = coalesce(qc_passed, FALSE),
    `QC status` = if_else(qc_passed, "QC passed", "Excluded"),
    `QC split reason` = case_when(
      qc_passed ~ "QC passed; included in Fig 2/3 and Fig 4 matched sample set",
      Species == "Neurospora crassa" ~ "Excluded from Fig 2/3 and Fig 4 species set",
      `Matched sample ID` == "AN_7" ~ "Excluded: AN_7 human-DNA contamination",
      TRUE ~ "Excluded: not in Fig 4 QC-passed/matched sample set"
    )
  ) %>%
  select(-qc_passed) %>%
  relocate(`QC status`, `QC split reason`, .after = `Matched sample ID`)

species_order <- c("Aspergillus niger", "Colletotrichum nagasakiense", "Neurospora crassa")
route_order <- c("Direct", "Indirect")

all_with_status <- all_with_status %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Route = factor(Route, levels = route_order)
  ) %>%
  arrange(Species, Route, `Matched sample ID`) %>%
  mutate(
    Species = as.character(Species),
    Route = as.character(Route)
  )

qc_passed <- all_with_status %>%
  filter(`QC status` == "QC passed")

excluded <- all_with_status %>%
  filter(`QC status` == "Excluded")

write_tsv(qc_passed, file.path(outdir, "QC_passed_samples.tsv"))
write_csv(qc_passed, file.path(outdir, "QC_passed_samples.csv"))
write_tsv(excluded, file.path(outdir, "Excluded_samples.tsv"))
write_csv(excluded, file.path(outdir, "Excluded_samples.csv"))
write_tsv(all_with_status, file.path(outdir, "All_samples_with_QC_status.tsv"))
write_csv(all_with_status, file.path(outdir, "All_samples_with_QC_status.csv"))

message("QC-passed samples:")
print(qc_passed %>% count(Species, Route, name = "n"), n = Inf)

message("Excluded samples:")
print(excluded %>% count(Species, Route, `QC split reason`, name = "n"), n = Inf)

message("Wrote outputs to: ", normalizePath(outdir))
