#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

data_root <- Sys.getenv(
  "PAPER_DATA_ROOT",
  unset = normalizePath(file.path(getwd(), "PaperDataFiles"), mustWork = FALSE)
)


repo_root <- file.path(data_root, "graphs", "standardisedpaperfigures")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) repo_root <- normalizePath(args[1], mustWork = TRUE)

outdir <- file.path(repo_root, "outputs", "fig2_fig3_memo_abcd")  # overridden below

paths <- list(
  direct = list(
    an_nn    = file.path(data_root, "assembly", "single_SAG_nonnorm", "AN_master_metrics.xlsx"),
    cntw_nn  = file.path(data_root, "assembly", "single_SAG_nonnorm", "CNTW_master_metrics.xlsx"),
    an_rm    = file.path(data_root, "reference_mapping", "refmap_metrics", "AN_shallow_direct_refmap_metrics_200k.csv"),
    cntw_rm  = file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_shallow_direct_refmap_metrics_200k.csv")
  ),
  indirect = list(
    an_nn    = file.path(data_root, "assembly", "single_SAG_nonnorm", "ANdeep_master_metrics.csv"),
    cntw_nn  = file.path(data_root, "assembly", "single_SAG_nonnorm", "CNTWdeep_master_metrics.csv"),
    an_rm    = file.path(data_root, "reference_mapping", "refmap_metrics", "AN_deep_indirect_refmap_metrics_200k.csv"),
    cntw_rm  = file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_deep_indirect_refmap_metrics_200k.csv")
  )
)

fig4_sample_paths <- list(
  direct = list(
    "A. niger" = file.path(data_root, "reference_mapping", "bin_qc", "AN_shallow_direct_bin_qc_10kb.csv"),
    "C. nagasakiense" = file.path(data_root, "reference_mapping", "bin_qc", "CNTW_shallow_direct_bin_qc_10kb.csv")
  ),
  indirect = list(
    "A. niger" = file.path(data_root, "reference_mapping", "bin_qc", "AN_deep_indirect_bin_qc_10kb.csv"),
    "C. nagasakiense" = file.path(data_root, "reference_mapping", "bin_qc", "CNTW_deep_indirect_bin_qc_10kb.csv")
  )
)

target_sample_counts <- list(
  direct = c("A. niger" = 175L, "C. nagasakiense" = 170L),
  indirect = c("A. niger" = 10L, "C. nagasakiense" = 15L)
)

sample_exclusions <- list(
  direct = list(
    "A. niger" = character(),
    "C. nagasakiense" = character()
  ),
  indirect = list(
    "A. niger" = character(),
    "C. nagasakiense" = character()
  )
)

cols <- list(
  direct = c("A. niger" = "#E1BE6A", "C. nagasakiense" = "#2D8B82"),
  indirect = c("A. niger" = "#1A85FF", "C. nagasakiense" = "#D41159")
)

display_group <- function(x) {
  recode(x, "A.niger" = "A. niger", "C.nagasakiense" = "C. nagasakiense", .default = x)
}

format_millions <- function(x) {
  paste0(sprintf("%.1f", x / 1e6), "M")
}

read_checked <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  read_csv(path, show_col_types = FALSE)
}

read_nn <- function(path, species_label) {
  ext <- tools::file_ext(path)
  df  <- if (ext == "xlsx") readxl::read_xlsx(path) else read_csv(path, show_col_types = FALSE)
  df %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    mutate(
      sample           = as.character(sample),
      group            = species_label,
      reads_after      = as.numeric(reads_after),
      busco_complete   = as.numeric(busco_complete),
      n50_kbp          = as.numeric(n50_bp) / 1000,
      total_length_mbp = as.numeric(totallength_bp) / 1e6,
      contig_count     = as.numeric(numcontigs)
    ) %>%
    select(sample, group, reads_after, busco_complete, n50_kbp, total_length_mbp, contig_count)
}

read_rm <- function(path, species_label) {
  read_csv(path, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    mutate(
      sample                 = as.character(sample),
      group                  = species_label,
      properly_paired_pct    = as.numeric(properly_paired_pct),
      mean_depth_covered_pos = as.numeric(mean_depth_all_pos),   # all-pos ref depth ~37x
      frac_covered_gt0       = as.numeric(frac_covered_gt0) * 100 # fraction -> %
    ) %>%
    select(sample, group, properly_paired_pct, mean_depth_covered_pos, frac_covered_gt0)
}

read_fig4_sample_set <- function(route, species_label) {
  path <- fig4_sample_paths[[route]][[species_label]]
  if (!file.exists(path)) stop("Missing Fig 4 sample file: ", path)

  samples <- read_csv(path, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    rename(sample = 1) %>%
    transmute(sample = as.character(sample)) %>%
    filter(!sample %in% sample_exclusions[[route]][[species_label]]) %>%
    distinct(sample) %>%
    pull(sample)

  expected_n <- target_sample_counts[[route]][[species_label]]
  if (length(samples) != expected_n) {
    stop(sprintf(
      "Fig 4 sample set mismatch for %s / %s: expected %d, found %d",
      route, species_label, expected_n, length(samples)
    ))
  }

  samples
}

load_route_data <- function(route) {
  fig4_samples <- bind_rows(lapply(names(cols[[route]]), function(species_label) {
    tibble(
      sample = read_fig4_sample_set(route, species_label),
      group = species_label
    )
  }))

  nn <- bind_rows(
    read_nn(paths[[route]]$an_nn,   names(cols[[route]])[1]),
    read_nn(paths[[route]]$cntw_nn, names(cols[[route]])[2])
  ) %>%
    filter(sample != "AN_7") %>%
    semi_join(fig4_samples, by = c("sample", "group"))

  rm <- bind_rows(
    read_rm(paths[[route]]$an_rm,   names(cols[[route]])[1]),
    read_rm(paths[[route]]$cntw_rm, names(cols[[route]])[2])
  ) %>%
    filter(sample != "AN_7") %>%
    semi_join(fig4_samples, by = c("sample", "group"))

  message(sprintf("[%s] nn: %d rows, samples: %s", route, nrow(nn), paste(head(nn$sample,3), collapse=", ")))
  message(sprintf("[%s] rm: %d rows, samples: %s", route, nrow(rm), paste(head(rm$sample,3), collapse=", ")))

  joined <- nn %>%
    inner_join(rm, by = c("sample", "group")) %>%
    mutate(group = factor(group, levels = names(cols[[route]])))

  joined_counts <- joined %>%
    distinct(sample, group) %>%
    count(group, name = "n")

  expected_counts <- tibble(
    group = factor(names(target_sample_counts[[route]]), levels = names(cols[[route]])),
    expected_n = as.integer(target_sample_counts[[route]])
  )

  count_check <- expected_counts %>%
    left_join(joined_counts, by = "group") %>%
    mutate(n = replace_na(n, 0L))

  if (any(count_check$n != count_check$expected_n)) {
    stop(
      "Joined Fig 2/3 sample counts do not match requested Fig 4 counts:\n",
      paste(
        sprintf("%s expected %d found %d", count_check$group, count_check$expected_n, count_check$n),
        collapse = "\n"
      )
    )
  }

  message(sprintf("[%s] joined: %d rows", route, nrow(joined)))
  joined
}

theme_pub <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1.4, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = base_size + 2, color = "black", hjust = 0.5, margin = margin(3, 0, 3, 0)),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      plot.margin = margin(7, 9, 9, 10)
    )
}

add_panel_tag <- function(p, tag, pos = c(0.035, 0.965)) {
  p +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(size = 13, face = "plain", color = "black"),
      plot.tag.position = pos
    )
}

row_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4.5, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.65),
      plot.margin = margin(2, 8, 2, 8)
    )
}

cor_label <- function(dat, x_col, y_col, log_x = FALSE) {
  dat %>%
    group_by(group) %>%
    summarise(
      r = suppressWarnings(cor(if (log_x) log10(.data[[x_col]]) else .data[[x_col]], .data[[y_col]], use = "complete.obs")),
      .groups = "drop"
    ) %>%
    mutate(txt = paste0(group, " r = ", sprintf("%.3f", r))) %>%
    pull(txt) %>%
    paste(collapse = "\n")
}

group_label_panel <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 1.12, label = label,
             size = 13 / .pt, fontface = "plain", color = "black",
             hjust = 0.5, vjust = 1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 4, 0))
}

make_metrics_panel <- function(dat, metric_spec, title, fill_cols, y_label = NULL, panel_tag = NULL, tag_pos = c(0.035, 0.965)) {
  panels <- lapply(seq_len(nrow(metric_spec)), function(i) {
    spec <- metric_spec[i, ]
    pd   <- tibble(group = dat$group, value = dat[[spec$column]] * spec$scale)
    y_br <- seq(0, spec$y_max, by = spec$y_step)

    pd$panel_label <- spec$label

    ggplot(pd, aes(x = group, y = value, fill = group)) +
      geom_boxplot(width = 0.58, outlier.shape = NA, color = "black", linewidth = 0.42, alpha = 0.72) +
      geom_jitter(shape = 21, width = 0.13, size = 1.25, alpha = 0.78, color = "black", stroke = 0.12) +
      scale_fill_manual(values = fill_cols, breaks = names(fill_cols)) +
      scale_x_discrete(labels = c("", "")) +
      scale_y_continuous(breaks = y_br, labels = label_comma(),
                         expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(ylim = c(0, max(y_br))) +
      facet_wrap(~ panel_label, nrow = 1) +
      labs(x = NULL, y = NULL) +
      theme_pub() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text   = element_text(size = 10.4, color = "black", hjust = 0.5,
                                    lineheight = 0.92, margin = margin(3, 0, 3, 0)),
        plot.margin  = margin(7, 4, 4, 4)
      )
  })

  if (!is.null(panel_tag)) panels[[1]] <- add_panel_tag(panels[[1]], panel_tag, pos = tag_pos)

  wrap_plots(panels, nrow = 1)
}

make_depth_coverage_panel <- function(dat, fill_cols) {
  label <- cor_label(dat, "mean_depth_covered_pos", "frac_covered_gt0")

  ggplot(dat, aes(x = mean_depth_covered_pos, y = frac_covered_gt0, fill = group)) +
    geom_point(shape = 21, size = 1.85, alpha = 0.84, color = "black", stroke = 0.15) +
    geom_smooth(aes(group = group), method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.65) +
    annotate("text", x = Inf, y = -Inf, label = label, hjust = 1.03, vjust = -0.4, size = 3.1, lineheight = 1.0) +
    facet_wrap(~ "Depth vs genome covered", nrow = 1) +
    scale_fill_manual(values = fill_cols, breaks = names(fill_cols)) +
    scale_x_continuous(breaks = pretty_breaks(n = 5), labels = label_comma(), expand = expansion(mult = c(0.04, 0.06))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Mean reference depth", y = "Genome covered (%)") +
    theme_pub() +
    theme(
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(7, 9, 4, 6)
    )
}

make_reads_busco_panel <- function(dat, fill_cols) {
  label <- cor_label(dat, "reads_after", "busco_complete", log_x = TRUE)
  x_breaks <- pretty(c(0, max(dat$reads_after, na.rm = TRUE)), n = 5)
  x_breaks <- x_breaks[x_breaks > 0]

  ggplot(dat, aes(x = reads_after, y = busco_complete, fill = group)) +
    geom_point(shape = 21, size = 1.85, alpha = 0.84, color = "black", stroke = 0.15) +
    geom_smooth(aes(group = group), method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.65) +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.05, vjust = 1.3, size = 3.1, lineheight = 1.0) +
    facet_wrap(~ "Reads vs BUSCO", nrow = 1) +
    scale_fill_manual(values = fill_cols, breaks = names(fill_cols)) +
    scale_x_continuous(breaks = x_breaks, labels = format_millions(x_breaks), expand = expansion(mult = c(0.04, 0.06))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Reads after trimming", y = "BUSCO completeness (%)") +
    theme_pub() +
    theme(
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(7, 9, 4, 6)
    )
}

make_reads_n50_panel <- function(dat, fill_cols) {
  label <- cor_label(dat, "reads_after", "n50_kbp", log_x = TRUE)
  x_breaks <- pretty(c(0, max(dat$reads_after, na.rm = TRUE)), n = 5)
  x_breaks <- x_breaks[x_breaks > 0]

  ggplot(dat, aes(x = reads_after, y = n50_kbp, fill = group)) +
    geom_point(shape = 21, size = 1.85, alpha = 0.84, color = "black", stroke = 0.15) +
    geom_smooth(aes(group = group), method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.65) +
    annotate("text", x = Inf, y = Inf, label = label, hjust = 1.03, vjust = 1.08, size = 3.1, lineheight = 1.0) +
    facet_wrap(~ "Reads vs N50", nrow = 1) +
    scale_fill_manual(values = fill_cols, breaks = names(fill_cols)) +
    scale_x_continuous(breaks = x_breaks, labels = format_millions(x_breaks), expand = expansion(mult = c(0.04, 0.06))) +
    scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(n = 5), expand = expansion(mult = c(0.03, 0.08))) +
    labs(x = "Reads after trimming", y = "N50 (kbp)") +
    theme_pub()
}

legend_plot <- function(dat, fill_cols) {
  counts <- dat %>%
    distinct(sample, group) %>%
    count(group, name = "n") %>%
    mutate(label = paste0(group, " (n=", n, ")"))

  legend_dat <- tibble(
    group = factor(names(fill_cols), levels = names(fill_cols)),
    x = c(0.20, 0.50)[seq_along(fill_cols)]
  ) %>%
    left_join(counts, by = "group")

  ggplot() +
    geom_rect(
      data = legend_dat,
      aes(xmin = x - 0.055, xmax = x - 0.015, ymin = 0.36, ymax = 0.64, fill = group),
      color = "black",
      linewidth = 0.25
    ) +
    geom_point(
      data = legend_dat,
      aes(x = x - 0.035, y = 0.5, fill = group),
      shape = 21,
      color = "black",
      size = 2.45,
      stroke = 0.2
    ) +
    geom_text(
      data = legend_dat,
      aes(x = x, y = 0.5, label = label),
      hjust = 0,
      size = 4.0,
      color = "black"
    ) +
    geom_segment(aes(x = 0.78, xend = 0.85, y = 0.5, yend = 0.5), color = "black", linetype = "dashed", linewidth = 0.7) +
    geom_text(aes(x = 0.87, y = 0.5, label = "Linear fit"), hjust = 0, size = 4.0, color = "black") +
    scale_fill_manual(values = fill_cols, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))
}

outdir <- file.path(repo_root, "outputs", "fig2_fig3_memo_abcd")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

build_figure <- function(route, file_prefix) {
  dat       <- load_route_data(route)
  fill_cols <- cols[[route]]

  # Panel A: read/reference-coverage metrics (4 sub-panels)
  read_metrics <- tibble::tribble(
    ~column,                   ~label,                ~scale, ~y_max, ~y_step,
    "reads_after",             "Total reads\n(M)",    1e-6,   15,     3,
    "properly_paired_pct",     "Mapped\nreads (%)",   1,      100,    20,
    "mean_depth_covered_pos",  "Mean depth\n(×)",     1,      50,     10,
    "frac_covered_gt0",        "Genome\ncovered (%)", 1,      100,    20
  )

  # Panel B: assembly metrics (3 sub-panels: length, contig count, N50)
  assembly_metrics <- tibble::tribble(
    ~column,            ~label,                   ~scale, ~y_max, ~y_step,
    "total_length_mbp", "Assembly\nlength (Mbp)",  1,      25,     5,
    "contig_count",     "Contig\nnumber",           1,      5000,   1000,
    "n50_kbp",          "N50 (kbp)\n",              1,      40,     8
  )

  p_a <- make_metrics_panel(dat, read_metrics,     "Read metrics",     fill_cols)
  p_b <- make_metrics_panel(dat, assembly_metrics, "Assembly metrics", fill_cols)
  p_c <- add_panel_tag(make_depth_coverage_panel(dat, fill_cols), "C")
  p_d <- add_panel_tag(make_reads_busco_panel(dat, fill_cols), "D") +
    theme(plot.margin = margin(7, 9, 4, 20))

  # Row 1: [A label] | A panels | spacer | [B label] | B panels
  row_1 <- wrap_plots(
    group_label_panel("A"), p_a, plot_spacer(), group_label_panel("B"), p_b,
    ncol = 5, widths = c(0.045, 1.48, 0.10, 0.045, 1.12)
  )

  # Row 2: C (depth vs covered) | D (reads vs BUSCO)
  row_2 <- p_c + p_d +
    plot_layout(ncol = 2, widths = c(1, 1))

  fig <- wrap_plots(
    list(row_1, row_2, legend_plot(dat, fill_cols)),
    ncol    = 1,
    heights = c(1.05, 1.0, 0.14)
  )

  ggsave(file.path(outdir, paste0(file_prefix, ".png")),
         fig, width = 12.0, height = 7.5, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(outdir, paste0(file_prefix, ".pdf")),
         fig, width = 12.0, height = 7.5, bg = "white", limitsize = FALSE)

  write_csv(dat, file.path(outdir, paste0(file_prefix, "_input_data.csv")))
  message("Saved: ", file_prefix)
}

build_figure("direct",   "fig2_direct_route_memo_abcd")
build_figure("indirect", "fig3_indirect_route_memo_abcd")

message("Done. Outputs written to: ", normalizePath(outdir))
