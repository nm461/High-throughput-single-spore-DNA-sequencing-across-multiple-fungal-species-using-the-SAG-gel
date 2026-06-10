#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
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

out_dir <- file.path(repo_root, "outputs", "UPDATED_suppS4_NC_memo_publication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nc_colour <- "#1A85FF"

read_csv_checked <- function(path) {
  if (!file.exists(path)) stop("Missing required data file: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

format_millions <- function(x) {
  paste0(sprintf("%.1f", x / 1e6), "M")
}

theme_pub <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x = element_text(size = base_size + 1, color = "black", margin = margin(t = 7)),
      axis.title.y = element_text(size = base_size + 1, color = "black", margin = margin(r = 8)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.5, face = "plain", color = "black", margin = margin(2, 0, 2, 0)),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.75),
      axis.ticks = element_line(color = "black", linewidth = 0.75),
      axis.ticks.length = unit(0.12, "cm"),
      plot.margin = margin(7, 9, 9, 9),
      plot.tag = element_text(size = base_size + 3, face = "plain", color = "black")
    )
}

add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(
      plot.tag.position = c(0.01, 1.04),
      plot.margin = margin(14, 9, 9, 14)
    )
}

panel_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.98, label = label, size = 5.0, color = "black", hjust = 0.5, vjust = 1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(4, 0, 0, 0))
}

tag_gutter <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.98, label = label, size = 5.0, color = "black", hjust = 0.5, vjust = 1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(7, 0, 0, 0))
}

metric_box <- function(dat, metric_order, title, panel_tag = NULL) {
  metric_order <- gsub("\\\\n", "\n", metric_order)
  dat <- dat %>%
    mutate(
      metric = gsub("\\\\n", "\n", metric),
      metric = factor(metric, levels = metric_order)
    )

  axis_specs <- list(
    "Total reads\n(M)" = list(max = 15, step = 3),
    "Mapped reads\n(%)" = list(max = 100, step = 20),
    "Mean depth\n(×)" = list(max = 250, step = 50),
    "Genome covered\n(%)" = list(max = 100, step = 20),
    "Assembly\nlength (Mbp)" = list(max = 25, step = 5),
    "Contig count" = list(max = 5000, step = 1000),
    "N50 (kbp)\n" = list(max = 40, step = 8)
  )

  panels <- lapply(seq_along(metric_order), function(i) {
    metric_name <- metric_order[[i]]
    pd <- dat %>% filter(metric == metric_name)
    spec <- axis_specs[[metric_name]]
    if (is.null(spec)) {
      ymax <- max(pd$value, 0, na.rm = TRUE)
      spec <- list(max = ymax * 1.1, step = ymax / 4)
    }
    breaks <- seq(0, spec$max, by = spec$step)

    p <- ggplot(pd, aes(x = "N. crassa", y = value)) +
      geom_boxplot(width = 0.44, fill = alpha(nc_colour, 0.45), color = "black", linewidth = 0.45, outlier.shape = NA) +
      geom_jitter(width = 0.09, size = 1.8, alpha = 0.85, color = nc_colour) +
      scale_y_continuous(breaks = breaks, labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(ylim = c(0, spec$max)) +
      labs(x = NULL, y = NULL, title = metric_name) +
      theme_pub(11.5) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12, face = "plain", color = "black", hjust = 0.5, lineheight = 0.92),
        plot.margin = margin(7, 6, 7, if (i == 1) 18 else 6)
      )

    if (i == 1 && !is.null(panel_tag)) {
      p <- p +
        labs(tag = panel_tag) +
        theme(
          plot.tag = element_text(size = 14.5, face = "plain", color = "black"),
          plot.tag.position = c(-0.03, 1.0)
        )
    }

    p
  })

  wrap_plots(panels, nrow = 1)
}

scatter_panel <- function(dat, x_col, y_col, title, x_lab, y_lab, tag, x_labels = NULL, y_limits = NULL) {
  r_val <- suppressWarnings(cor(dat[[x_col]], dat[[y_col]], use = "complete.obs"))
  p <- ggplot(dat, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(shape = 21, size = 2.0, fill = nc_colour, color = "black", stroke = 0.18, alpha = 0.88) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = Inf, y = Inf, label = paste0("N. crassa r = ", sprintf("%.3f", r_val)), hjust = 1.04, vjust = 1.12, size = 3.5, color = "black") +
    facet_wrap(~ factor(title, levels = title), nrow = 1) +
    labs(x = x_lab, y = y_lab) +
    theme_pub(11.5)
  if (!is.null(x_labels)) {
    p <- p + scale_x_continuous(labels = x_labels, breaks = pretty_breaks(n = 4))
  } else {
    p <- p + scale_x_continuous(breaks = pretty_breaks(n = 4))
  }
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 25), expand = expansion(mult = c(0, 0.05)))
  } else {
    p <- p + scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(n = 5), expand = expansion(mult = c(0, 0.05)))
  }
  if (is.null(tag) || tag == "") {
    p
  } else {
    add_tag(p, tag)
  }
}

readcov <- read_csv_checked(file.path(repo_root, "output", "suppS_NC_D", "panelG_combined_input_data.csv")) %>%
  mutate(
    sample = as.character(sample),
    properly_paired_pct = as.numeric(properly_paired_pct),
    frac_covered_gt0 = as.numeric(frac_covered_gt0),
    mean_depth_covered_pos = as.numeric(mean_depth_covered_pos)
  )

scat <- read_csv_checked(file.path(repo_root, "output", "suppS_NC_A", "panelA_combined_input_data.csv")) %>%
  mutate(
    sample = as.character(sample),
    reads_after = as.numeric(reads_after),
    busco_complete = as.numeric(busco_complete),
    n50_bp = as.numeric(n50_bp),
    n50_kbp = n50_bp / 1000
  )

asm <- read_csv_checked(file.path(repo_root, "output", "suppS_NC_BC", "A3_input_data_used.csv")) %>%
  mutate(
    sample = as.character(sample),
    total_length_mbp = as.numeric(total_length_mbp),
    contig_count = as.numeric(contig_count),
    n50_kbp = as.numeric(n50_kbp)
  )

read_panel_dat <- readcov %>%
  left_join(scat %>% select(sample, reads_after), by = "sample") %>%
  pivot_longer(
    cols = c(reads_after, properly_paired_pct, mean_depth_covered_pos, frac_covered_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(metric == "reads_after", value / 1e6, value),
    metric = recode(
      metric,
      reads_after = "Total reads\\n(M)",
      properly_paired_pct = "Mapped reads\\n(%)",
      mean_depth_covered_pos = "Mean depth\n(×)",
      frac_covered_gt0 = "Genome covered\\n(%)"
    ),
    panel = "Read and reference-coverage metrics"
  ) %>%
  filter(is.finite(value))

assembly_panel_dat <- asm %>%
  pivot_longer(
    cols = c(total_length_mbp, contig_count, n50_kbp),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      total_length_mbp = "Assembly\\nlength (Mbp)",
      contig_count = "Contig count",
      n50_kbp = "N50 (kbp)\\n"
    ),
    panel = "Assembly metrics"
  ) %>%
  filter(is.finite(value))

p_a <- (tag_gutter("A") | metric_box(
    read_panel_dat,
    c("Total reads\\n(M)", "Mapped reads\\n(%)", "Mean depth\\n(×)", "Genome covered\\n(%)"),
    "Read and reference-coverage metrics"
  )) + plot_layout(widths = c(0.055, 1))

p_b <- (tag_gutter("B") | metric_box(
    assembly_panel_dat,
    c("Assembly\\nlength (Mbp)", "Contig count", "N50 (kbp)\\n"),
    "Assembly metrics"
  )) + plot_layout(widths = c(0.055, 1))

p_c <- scatter_panel(
  readcov,
  "mean_depth_covered_pos",
  "frac_covered_gt0",
  "Mean depth vs genome\ncovered percentage",
  "Estimated mean reference depth",
  "Genome covered (%)",
  "",
  y_limits = c(0, 100)
)

p_d <- scatter_panel(
  scat,
  "reads_after",
  "busco_complete",
  "Reads vs BUSCO",
  "Reads after\ntrimming",
  "BUSCO\ncompleteness (%)",
  "D",
  x_labels = format_millions,
  y_limits = c(0, 100)
)

gini_dat <- read_csv_checked(file.path(repo_root, "output", "suppS_NC_gini", "gini_by_group_data.csv")) %>%
  mutate(gini = as.numeric(gini), panel = "Gini coefficient distribution")

gini_summary <- gini_dat %>%
  summarise(n = n(), median_gini = median(gini, na.rm = TRUE), .groups = "drop")

p_e_base <- ggplot(gini_dat, aes(x = "N. crassa", y = gini)) +
  geom_jitter(width = 0.10, size = 2.0, alpha = 0.88, color = nc_colour) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, color = "black", linewidth = 0.5) +
  annotate("label", x = 1.18, y = 0.48, label = sprintf("Median Gini\nN. crassa  %.2f", gini_summary$median_gini), hjust = 0, size = 2.75, fill = "white", color = "black") +
  facet_wrap(~ panel, nrow = 1) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Sample\ngroup", y = "Gini coefficient") +
  theme_pub(11.5)
p_e <- p_e_base

lorenz <- read_csv_checked(file.path(repo_root, "outputs", "suppNC_fig4AB_gini_lorenz_publication", "suppNC_fig4B_lorenz_summary.csv"))

p_f_base <- ggplot(lorenz, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.65) +
  geom_line(color = nc_colour, linewidth = 1.0) +
  facet_wrap(~ factor("Lorenz curve", levels = "Lorenz curve"), nrow = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  labs(x = "Cumulative fraction of genome bins", y = "Cumulative fraction\nof total coverage") +
  theme_pub(11.5)
p_f <- add_tag(p_f_base, "F")

coassembly <- read_csv_checked(file.path(repo_root, "outputs", "suppNC_coassembly_publication", "suppNC_coassembly_publication_plot_data.csv")) %>%
  filter(metric == "Completeness") %>%
  mutate(sag = as.numeric(sag), value = as.numeric(value))

coassembly_label <- coassembly %>%
  filter(sag == max(sag, na.rm = TRUE)) %>%
  mutate(label = sprintf("%.1f%%", value), label_x = sag + 0.12, label_y = value + 1.2)

p_g_base <- ggplot(coassembly, aes(x = sag, y = value)) +
  geom_line(color = nc_colour, linewidth = 1.0) +
  geom_point(color = nc_colour, size = 2.0) +
  geom_text(data = coassembly_label, aes(x = label_x, y = label_y, label = label), hjust = 0, color = nc_colour, size = 3.2) +
  facet_wrap(~ factor("Co-assembly BUSCO completeness", levels = "Co-assembly BUSCO completeness"), nrow = 1) +
  scale_x_continuous(breaks = 2:12, limits = c(1.7, 13.3)) +
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of SAGs co-assembled", y = "BUSCO completeness") +
  theme_pub(11.5)
p_g <- p_g_base

legend_plot <- ggplot() +
  geom_segment(aes(x = 0.30, xend = 0.38, y = 0.5, yend = 0.5), color = nc_colour, linewidth = 1.1, lineend = "round") +
  geom_point(aes(x = 0.34, y = 0.5), color = nc_colour, size = 2.6) +
  geom_text(aes(x = 0.40, y = 0.5, label = sprintf("N. crassa (n=%s)", gini_summary$n)), hjust = 0, size = 4.0, color = "black") +
  geom_segment(aes(x = 0.62, xend = 0.70, y = 0.5, yend = 0.5), color = "black", linetype = "dashed", linewidth = 0.75) +
  geom_text(aes(x = 0.72, y = 0.5, label = "Linear fit / equality line"), hjust = 0, size = 4.0, color = "black") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

row_cd <- (tag_gutter("C") | (p_c | p_d)) +
  plot_layout(widths = c(0.055, 1))
row_ef <- (tag_gutter("E") | (p_e | p_f)) +
  plot_layout(widths = c(0.055, 1))
row_g <- (tag_gutter("G") | p_g) +
  plot_layout(widths = c(0.055, 1))

fig <- p_a / p_b / row_cd / row_ef / row_g / legend_plot +
  plot_layout(heights = c(0.75, 0.75, 1, 1, 0.95, 0.12))

ggsave(
  file.path(out_dir, "UPDATED_suppS4_NC_memo_publication.png"),
  fig,
  width = 8.8,
  height = 15.0,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(out_dir, "UPDATED_suppS4_NC_memo_publication.pdf"),
  fig,
  width = 8.8,
  height = 15.0,
  bg = "white",
  limitsize = FALSE
)

write_csv(read_panel_dat, file.path(out_dir, "UPDATED_suppS4_NC_read_reference_metrics_data.csv"))
write_csv(assembly_panel_dat, file.path(out_dir, "UPDATED_suppS4_NC_assembly_metrics_data.csv"))
write_csv(coassembly, file.path(out_dir, "UPDATED_suppS4_NC_coassembly_completeness_data.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir))
