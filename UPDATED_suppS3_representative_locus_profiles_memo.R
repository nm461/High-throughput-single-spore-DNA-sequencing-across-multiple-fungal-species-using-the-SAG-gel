#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(grid)
  library(scales)
})

data_root <- Sys.getenv(
  "PAPER_DATA_ROOT",
  unset = normalizePath(file.path(getwd(), "PaperDataFiles"), mustWork = FALSE)
)


repo_root <- file.path(data_root, "graphs", "standardisedpaperfigures")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) repo_root <- normalizePath(args[1], mustWork = TRUE)
if (length(args) >= 2) data_root <- normalizePath(args[2], mustWork = TRUE)

cluster_root <- file.path(data_root, "cluster_inputs")
outdir <- file.path(repo_root, "outputs", "UPDATED_suppS3_representative_locus_profiles_memo")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

coverage_threshold <- 1.0
bin_file <- "bin_qc_10kb.tsv.gz"

cols <- c(
  "A. niger|Direct route" = "#E1BE6A",
  "A. niger|Indirect route" = "#1A85FF",
  "C. nagasakiense|Direct route" = "#2D8B82",
  "C. nagasakiense|Indirect route" = "#D41159"
)

included_loci <- tibble::tribble(
  ~figure_id, ~row_group, ~panel, ~species, ~locus_label, ~contig, ~feature_start, ~feature_end, ~plot_start, ~plot_end,
  "suppS3", "A     C. nagasakiense BUSCO loci", "i", "C. nagasakiense", "BUSCO 244066at4751", "cntw_is_4", 2200778L, 2202195L, 2170778L, 2232195L,
  "suppS3", "A     C. nagasakiense BUSCO loci", "ii", "C. nagasakiense", "BUSCO 445210at4751", "cntw_is_10", 419915L, 420155L, 389915L, 450155L,
  "suppS3", "B     C. nagasakiense coverage-anomalous loci", "i", "C. nagasakiense", "cntw_is_7:2610-2620 kb", "cntw_is_7", 2610000L, 2620000L, 2580000L, 2650000L,
  "suppS3", "B     C. nagasakiense coverage-anomalous loci", "ii", "C. nagasakiense", "cntw_is_7:2490-2500 kb", "cntw_is_7", 2490000L, 2500000L, 2460000L, 2530000L,
  "suppS3", "C     A. niger BUSCO loci", "i", "A. niger", "BUSCO 471205at4751", "NT_166519.1", 258307L, 258635L, 228307L, 288635L,
  "suppS3", "C     A. niger BUSCO loci", "ii", "A. niger", "BUSCO 122824at4751", "NT_166519.1", 3476809L, 3478371L, 3446809L, 3508371L,
  "suppS3", "D     A. niger coverage-anomalous loci", "i", "A. niger", "NT_166519.1:3410-3420 kb", "NT_166519.1", 3410000L, 3420000L, 3380000L, 3450000L,
  "suppS3", "D     A. niger coverage-anomalous loci", "ii", "A. niger", "NT_166526.1:860-870 kb", "NT_166526.1", 860000L, 870000L, 830000L, 900000L
)

included_loci <- included_loci %>%
  mutate(
    locus_key = paste(figure_id, row_group, species, locus_label, sep = "|"),
    title = paste0(panel, "     ", locus_label)
  )

promoted_to_fig4c <- c(
  "NT_166526.1:860-870 kb",
  "cntw_is_7:2610-2620 kb"
)

included_loci <- included_loci %>%
  filter(!locus_label %in% promoted_to_fig4c)

read_bins_for_method <- function(root, method_label, species_label, loci) {
  sample_files <- list.files(root, pattern = bin_file, recursive = TRUE, full.names = TRUE)
  if (length(sample_files) == 0) {
    stop("No ", bin_file, " files found in: ", root)
  }

  pieces <- vector("list", length(sample_files))
  for (i in seq_along(sample_files)) {
    fp <- sample_files[[i]]
    sample <- basename(dirname(fp))
    df <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)

    needed <- c("contig", "start", "end", "mean_depth", "gc_frac")
    missing <- setdiff(needed, names(df))
    if (length(missing) > 0) {
      stop("Missing columns in ", fp, ": ", paste(missing, collapse = ", "))
    }

    df <- df %>%
      mutate(
        start = as.integer(start),
        end = as.integer(end),
        mean_depth = as.numeric(mean_depth),
        gc_frac = as.numeric(gc_frac)
      )

    denom <- mean(df$mean_depth[is.finite(df$mean_depth)], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1

    locus_rows <- lapply(seq_len(nrow(loci)), function(j) {
      loc <- loci[j, ]
      df %>%
        filter(
          contig == loc$contig,
          end > loc$plot_start,
          start < loc$plot_end
        ) %>%
        transmute(
          sample = sample,
          species = species_label,
          method = method_label,
          figure_id = loc$figure_id,
          row_group = loc$row_group,
          locus_key = loc$locus_key,
          locus_label = loc$locus_label,
          panel = loc$panel,
          title = loc$title,
          contig = contig,
          start = start,
          end = end,
          mid_kb = ((start + end) / 2) / 1000,
          feature_start_kb = loc$feature_start / 1000,
          feature_end_kb = loc$feature_end / 1000,
          plot_start_kb = loc$plot_start / 1000,
          plot_end_kb = loc$plot_end / 1000,
          gc_frac = gc_frac,
          norm_depth = mean_depth / denom,
          covered = as.numeric(mean_depth >= coverage_threshold)
        )
    })
    pieces[[i]] <- bind_rows(locus_rows)
  }

  bind_rows(pieces)
}

an_loci <- included_loci %>% filter(species == "A. niger")
cntw_loci <- included_loci %>% filter(species == "C. nagasakiense")

raw_dat <- bind_rows(
  read_bins_for_method(file.path(cluster_root, "AN_shallow"), "Direct route", "A. niger", an_loci),
  read_bins_for_method(file.path(cluster_root, "AN_deep"), "Indirect route", "A. niger", an_loci),
  read_bins_for_method(file.path(cluster_root, "CNTW_shallow"), "Direct route", "C. nagasakiense", cntw_loci),
  read_bins_for_method(file.path(cluster_root, "CNTW_deep"), "Indirect route", "C. nagasakiense", cntw_loci)
)

plot_dat <- raw_dat %>%
  group_by(
    figure_id, row_group, species, method, locus_key, locus_label, panel, title, contig,
    start, end, mid_kb, feature_start_kb, feature_end_kb, plot_start_kb, plot_end_kb
  ) %>%
  summarise(
    gc_frac = first(na.omit(gc_frac)),
    mean_norm_depth = mean(norm_depth, na.rm = TRUE),
    prop_covered = mean(covered, na.rm = TRUE),
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  mutate(
    colour_key = paste(species, method, sep = "|"),
    title = factor(title, levels = included_loci$title)
  )

sample_counts <- plot_dat %>%
  distinct(species, method, n_samples, colour_key) %>%
  arrange(match(colour_key, names(cols)))

theme_locus <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.55),
      axis.ticks = element_line(color = "black", linewidth = 0.55),
      plot.title = element_text(size = base_size + 2, face = "plain", color = "black", hjust = 0.02),
      plot.margin = margin(2, 5, 2, 5)
    )
}

inner_x_breaks <- function(lims) {
  br <- pretty(lims, n = 4)
  br <- br[br > lims[1] & br < lims[2]]
  if (length(br) < 2) br <- pretty(lims, n = 3)
  br
}

make_locus_panel <- function(loc) {
  dat <- plot_dat %>% filter(locus_key == loc$locus_key)
  x_lims <- c(loc$plot_start / 1000, loc$plot_end / 1000)
  dat <- dat %>% filter(mid_kb >= x_lims[1], mid_kb <= x_lims[2])
  highlight_base <- data.frame(
    xmin = loc$feature_start / 1000,
    xmax = loc$feature_end / 1000
  )
  gc_dat <- dat %>%
    distinct(mid_kb, gc_frac, .keep_all = TRUE) %>%
    arrange(mid_kb)
  y_max <- max(dat$mean_norm_depth, na.rm = TRUE)
  depth_lim_hi <- max(1.5, y_max * 1.10)
  gc_lo <- min(gc_dat$gc_frac, na.rm = TRUE)
  gc_hi <- max(gc_dat$gc_frac, na.rm = TRUE)
  gc_pad <- max(0.015, (gc_hi - gc_lo) * 0.25)
  gc_lims <- c(max(0, gc_lo - gc_pad), min(1, gc_hi + gc_pad))
  highlight_gc <- mutate(highlight_base, ymin = gc_lims[1], ymax = gc_lims[2])
  highlight_depth <- mutate(highlight_base, ymin = 0, ymax = depth_lim_hi)
  highlight_cov <- mutate(highlight_base, ymin = 0, ymax = 1.02)

  p_gc <- ggplot(gc_dat, aes(x = mid_kb, y = gc_frac)) +
    geom_rect(data = highlight_gc, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey90", alpha = 0.35) +
    geom_line(color = "grey35", linewidth = 0.55) +
    scale_x_continuous(limits = x_lims, breaks = inner_x_breaks(x_lims), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(limits = gc_lims, breaks = pretty_breaks(n = 2), expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = "GC", title = loc$title) +
    theme_locus() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_depth <- ggplot(dat, aes(x = mid_kb, y = mean_norm_depth, color = colour_key, fill = colour_key)) +
    geom_rect(data = highlight_depth, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey90", alpha = 0.35) +
    geom_area(alpha = 0.10, linewidth = 0, position = "identity") +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, color = "grey60", linetype = "dotted", linewidth = 0.35) +
    annotate("text", x = x_lims[2], y = y_max, label = sprintf("[0-%.0fx]", y_max), hjust = 1, vjust = 1, size = 2.4, color = "grey30") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(limits = x_lims, breaks = inner_x_breaks(x_lims), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(limits = c(0, depth_lim_hi), breaks = pretty_breaks(n = 3), expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Norm. depth") +
    guides(color = "none", fill = "none") +
    theme_locus() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_cov <- ggplot(dat, aes(x = mid_kb, y = prop_covered, color = colour_key, fill = colour_key)) +
    geom_rect(data = highlight_cov, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey90", alpha = 0.35) +
    geom_area(alpha = 0.10, linewidth = 0, position = "identity") +
    geom_line(linewidth = 0.65) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(limits = x_lims, breaks = inner_x_breaks(x_lims), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Position on reference (kb)", y = "Cov.") +
    guides(color = "none", fill = "none") +
    theme_locus()

  p_gc / p_depth / p_cov + plot_layout(heights = c(0.82, 2.35, 0.62))
}

legend_plot <- function() {
  legend_dat <- sample_counts %>%
    mutate(
      label = sprintf("%s: %s (n=%s)", method, species, n_samples),
      x = c(0.08, 0.55, 0.08, 0.55)[seq_len(n())],
      y = c(0.68, 0.68, 0.28, 0.28)[seq_len(n())]
    )

  ggplot(legend_dat) +
    geom_segment(aes(x = x, xend = x + 0.08, y = y, yend = y, color = colour_key), linewidth = 0.9, lineend = "round") +
    geom_text(aes(x = x + 0.10, y = y, label = label), hjust = 0, vjust = 0.5, size = 3.4, color = "black") +
    scale_color_manual(values = cols, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 5, 2, 5))
}

row_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4.8, fontface = "plain", color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
      plot.margin = margin(3, 5, 3, 5)
    )
}

row_title <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4.8, fontface = "plain", color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 5, 2, 5))
}

panels_for <- function(figure, group) {
  locs <- included_loci %>%
    filter(figure_id == figure, row_group == group)
  split(locs, seq_len(nrow(locs))) %>%
    lapply(make_locus_panel)
}

make_two_row_figure <- function(figure, top_group, bottom_group, ncol, legend = TRUE) {
  parts <- row_strip(top_group) /
    plot_spacer() /
    wrap_plots(panels_for(figure, top_group), ncol = ncol) /
    plot_spacer() /
    row_strip(bottom_group) /
    plot_spacer() /
    wrap_plots(panels_for(figure, bottom_group), ncol = ncol)

  if (legend) {
    parts <- parts / legend_plot() +
      plot_layout(heights = c(0.08, 0.035, 1, 0.12, 0.08, 0.035, 1, 0.14))
  } else {
    parts <- parts + plot_layout(heights = c(0.08, 0.035, 1, 0.12, 0.08, 0.035, 1))
  }
  parts
}

make_suppS3_figure <- function() {
  busco_loci <- included_loci %>%
    filter(grepl("BUSCO", row_group)) %>%
    arrange(species, panel)
  anomalous_loci <- included_loci %>%
    filter(grepl("coverage-anomalous", row_group)) %>%
    arrange(species, panel)

  busco_panels <- split(busco_loci, seq_len(nrow(busco_loci))) %>%
    lapply(make_locus_panel)
  anomalous_panels <- split(anomalous_loci, seq_len(nrow(anomalous_loci))) %>%
    lapply(make_locus_panel)

  row_title("A. niger BUSCO loci") /
    (busco_panels[[1]] | busco_panels[[2]]) /
    patchwork::plot_spacer() /
    row_title("C. nagasakiense BUSCO loci") /
    (busco_panels[[3]] | busco_panels[[4]]) /
    patchwork::plot_spacer() /
    row_title("Coverage-anomalous loci") /
    wrap_plots(anomalous_panels, ncol = 2) /
    legend_plot() +
    plot_layout(heights = c(0.10, 1, 0.08, 0.10, 1, 0.18, 0.10, 1, 0.16))
}

fig_s3 <- make_suppS3_figure()

ggsave(
  file.path(outdir, "UPDATED_suppS3_representative_locus_profiles_memo.png"),
  fig_s3,
  width = 9.8,
  height = 15.2,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(outdir, "UPDATED_suppS3_representative_locus_profiles_memo.pdf"),
  fig_s3,
  width = 9.8,
  height = 15.2,
  bg = "white",
  limitsize = FALSE
)

write_csv(included_loci, file.path(outdir, "supp_locus_profile_variants_included.csv"))
write_csv(plot_dat, file.path(outdir, "supp_locus_profile_variants_plot_data.csv"))

message("Done. Outputs written to: ", normalizePath(outdir))
