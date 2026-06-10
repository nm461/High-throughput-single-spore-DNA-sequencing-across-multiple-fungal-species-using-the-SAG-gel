#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tools)
  library(grid)
  library(patchwork)
  library(scales)
})

data_root <- Sys.getenv(
  "PAPER_DATA_ROOT",
  unset = normalizePath(file.path(getwd(), "PaperDataFiles"), mustWork = FALSE)
)


has_readxl <- requireNamespace("readxl", quietly = TRUE)

file_an_direct_master       <- file.path(data_root, "reference_mapping", "refmap_metrics", "AN_shallow_direct_refmap_metrics_200k.csv")
file_an_direct_bin          <- file.path(data_root, "reference_mapping", "bin_qc", "AN_shallow_direct_bin_qc_10kb.csv")
file_cntw_direct_master     <- file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_shallow_direct_refmap_metrics_200k.csv")
file_cntw_direct_bin        <- file.path(data_root, "reference_mapping", "bin_qc", "CNTW_shallow_direct_bin_qc_10kb.csv")
file_an_indirect_master     <- file.path(data_root, "reference_mapping", "refmap_metrics", "AN_deep_indirect_refmap_metrics_200k.csv")
file_an_indirect_bin        <- file.path(data_root, "reference_mapping", "bin_qc", "AN_deep_indirect_bin_qc_10kb.csv")
file_cntw_indirect_master   <- file.path(data_root, "reference_mapping", "refmap_metrics", "CNTW_deep_indirect_refmap_metrics_200k.csv")
file_cntw_indirect_bin      <- file.path(data_root, "reference_mapping", "bin_qc", "CNTW_deep_indirect_bin_qc_10kb.csv")
file_an_direct_lorenz       <- file.path(data_root, "reference_mapping", "lorenz", "AN_shallow_direct_lorenz_curve_10kb.csv")
file_cntw_direct_lorenz     <- file.path(data_root, "reference_mapping", "lorenz", "CNTW_shallow_direct_lorenz_curve_10kb.csv")
file_an_indirect_lorenz     <- file.path(data_root, "reference_mapping", "lorenz", "AN_deep_indirect_lorenz_curve_10kb.csv")
file_cntw_indirect_lorenz   <- file.path(data_root, "reference_mapping", "lorenz", "CNTW_deep_indirect_lorenz_curve_10kb.csv")
outdir                      <- file.path(data_root, "graphs", "standardisedpaperfigures", "outputs", "UPDATED_fig4AB_gini_lorenz_memo_publication")
cluster_root                <- file.path(data_root, "cluster_inputs")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_any <- function(path) {
  ext <- tolower(file_ext(path))

  if (ext == "csv") {
    return(read_csv(path, show_col_types = FALSE))
  }

  if (ext %in% c("tsv", "txt")) {
    return(read_tsv(path, show_col_types = FALSE))
  }

  if (ext %in% c("xlsx", "xls")) {
    if (!has_readxl) {
      stop("Package 'readxl' is required for Excel files. Install with install.packages('readxl').")
    }
    return(readxl::read_excel(path))
  }

  stop(paste("Unsupported file type:", path))
}

clean_names_simple <- function(x) {
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  tolower(x)
}

check_required_cols <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required columns in ", label, ": ",
        paste(missing, collapse = ", "),
        "\nAvailable columns: ",
        paste(names(df), collapse = ", ")
      )
    )
  }
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
}

read_binqc_file <- function(path, group_label) {
  df <- read_any(path)
  names(df) <- clean_names_simple(names(df))

  required_cols <- c("sample", "gini_nonzero")
  check_required_cols(df, required_cols, basename(path))

  df %>%
    transmute(
      sample = as.character(sample),
      group = group_label,
      gini = suppressWarnings(as.numeric(gini_nonzero))
    ) %>%
    filter(!is.na(gini))
}

read_lorenz_file <- function(path, group_label) {
  df <- read_any(path)
  names(df) <- clean_names_simple(names(df))

  required_cols <- c(
    "sample",
    "binrank",
    "cumulativebinfraction",
    "cumulativedepthfraction"
  )
  check_required_cols(df, required_cols, basename(path))

  df %>%
    transmute(
      sample = as.character(sample),
      group = group_label,
      bin_rank = as.numeric(binrank),
      x = as.numeric(cumulativebinfraction),
      y = as.numeric(cumulativedepthfraction)
    ) %>%
    filter(!is.na(bin_rank), !is.na(x), !is.na(y))
}

invisible(lapply(
  c(
    file_an_direct_master, file_an_direct_bin,
    file_cntw_direct_master, file_cntw_direct_bin,
    file_an_indirect_master, file_an_indirect_bin,
    file_cntw_indirect_master, file_cntw_indirect_bin,
    file_an_direct_lorenz, file_cntw_direct_lorenz,
    file_an_indirect_lorenz, file_cntw_indirect_lorenz
  ),
  check_file_exists
))

group_order <- c("AN_direct", "CNTW_direct", "AN_indirect", "CNTW_indirect")

group_cols <- c(
  "AN_direct"     = "#E1BE6A",
  "CNTW_direct"   = "#2D8B82",
  "AN_indirect"   = "#1A85FF",
  "CNTW_indirect" = "#D41159"
)

gini_dat <- bind_rows(
  read_binqc_file(file_an_direct_bin,     "AN_direct"),
  read_binqc_file(file_cntw_direct_bin,   "CNTW_direct"),
  read_binqc_file(file_an_indirect_bin,   "AN_indirect"),
  read_binqc_file(file_cntw_indirect_bin, "CNTW_indirect")
) %>%
  mutate(
    group = factor(group, levels = group_order),
    panel = "Gini coefficient distributions"
  )

gini_group_summary <- gini_dat %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_gini = median(gini, na.rm = TRUE),
    .groups = "drop"
  )

n_lookup <- setNames(gini_group_summary$n, as.character(gini_group_summary$group))
gini_group_labels <- c(
  "AN_direct"     = sprintf("A. niger\nDirect route\n(n=%s)", n_lookup[["AN_direct"]]),
  "CNTW_direct"   = sprintf("C. nagasakiense\nDirect route\n(n=%s)", n_lookup[["CNTW_direct"]]),
  "AN_indirect"   = sprintf("A. niger\nIndirect route\n(n=%s)", n_lookup[["AN_indirect"]]),
  "CNTW_indirect" = sprintf("C. nagasakiense\nIndirect route\n(n=%s)", n_lookup[["CNTW_indirect"]])
)

p_to_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("***  p < 0.001")
  if (p < 0.01) return("**  p < 0.01")
  if (p < 0.05) return("*  p < 0.05")
  paste0("p = ", signif(p, 2))
}

p_an <- wilcox.test(
  gini ~ group,
  data = filter(gini_dat, group %in% c("AN_direct", "AN_indirect")),
  exact = FALSE
)$p.value

p_cntw <- wilcox.test(
  gini ~ group,
  data = filter(gini_dat, group %in% c("CNTW_direct", "CNTW_indirect")),
  exact = FALSE
)$p.value

bracket_data <- tibble::tibble(
  x = c(1, 2),
  xend = c(3, 4),
  y = c(0.925, 0.99),
  label = c(p_to_label(p_an), p_to_label(p_cntw))
)

median_text <- paste(
  "Median Gini",
  sprintf("A. niger Direct route   %.2f", gini_group_summary$median_gini[gini_group_summary$group == "AN_direct"]),
  sprintf("C. nagasakiense Direct route   %.2f", gini_group_summary$median_gini[gini_group_summary$group == "CNTW_direct"]),
  sprintf("A. niger Indirect route %.2f", gini_group_summary$median_gini[gini_group_summary$group == "AN_indirect"]),
  sprintf("C. nagasakiense Indirect route %.2f", gini_group_summary$median_gini[gini_group_summary$group == "CNTW_indirect"]),
  sep = "\n"
)

lorenz_dat <- bind_rows(
  read_lorenz_file(file_an_direct_lorenz,     "AN_direct"),
  read_lorenz_file(file_cntw_direct_lorenz,   "CNTW_direct"),
  read_lorenz_file(file_an_indirect_lorenz,   "AN_indirect"),
  read_lorenz_file(file_cntw_indirect_lorenz, "CNTW_indirect")
)

lorenz_summary <- lorenz_dat %>%
  group_by(group, bin_rank) %>%
  summarise(
    x = mean(x, na.rm = TRUE),
    y = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = group_order),
    panel = "Lorenz curves"
  )

lorenz_group_summary <- lorenz_dat %>%
  group_by(group) %>%
  summarise(
    n = n_distinct(sample),
    .groups = "drop"
  )

lorenz_group_labels <- c(
  "AN_direct"     = sprintf("Direct route: A. niger (n=%s)", lorenz_group_summary$n[lorenz_group_summary$group == "AN_direct"]),
  "CNTW_direct"   = sprintf("Direct route: C. nagasakiense (n=%s)", lorenz_group_summary$n[lorenz_group_summary$group == "CNTW_direct"]),
  "AN_indirect"   = sprintf("Indirect route: A. niger (n=%s)", lorenz_group_summary$n[lorenz_group_summary$group == "AN_indirect"]),
  "CNTW_indirect" = sprintf("Indirect route: C. nagasakiense (n=%s)", lorenz_group_summary$n[lorenz_group_summary$group == "CNTW_indirect"])
)

theme_paper <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title.x = element_text(
        size = 12, face = "plain", color = "black",
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 12, face = "plain", color = "black",
        margin = margin(r = 9)
      ),
      axis.text = element_text(size = 12, face = "plain", color = "black"),
      axis.text.x = element_text(size = 8.4, face = "plain", color = "black"),
      strip.text = element_text(
        size = 12,
        face = "plain",
        color = "black",
        margin = margin(2, 0, 2, 0)
      ),
      strip.background = element_rect(
        fill = NA,
        color = NA
      ),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "plain", color = "black"),
      legend.key = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.14, "cm"),
      plot.margin = margin(t = 10, r = 12, b = 16, l = 18)
    )
}

p_gini <- ggplot(gini_dat, aes(x = group, y = gini, color = group)) +
  geom_jitter(
    width = 0.13,
    size = 2.0,
    alpha = 0.85
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x, xend = x, y = y, yend = y - 0.025),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = xend, xend = xend, y = y, yend = y - 0.025),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    data = bracket_data,
    aes(x = (x + xend) / 2, y = y + 0.018, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 3.2
  ) +
  annotate(
    "label",
    x = 2.78,
    y = 0.34,
    label = median_text,
    hjust = 0,
    vjust = 1,
    size = 2.25,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(values = group_cols) +
  scale_x_discrete(labels = gini_group_labels) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Route and species",
    y = "Gini coefficient",
    tag = "A"
  ) +
  coord_cartesian(clip = "off") +
  theme_paper() +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "plain", color = "black"),
    plot.tag.position = c(-0.05, 1.0)
  )

p_lorenz <- ggplot(
  lorenz_summary,
  aes(x = x, y = y, color = group)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.7
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(
    values = group_cols,
    breaks = group_order,
    labels = lorenz_group_labels,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Cumulative fraction of genome bins",
    y = "Cumulative fraction of total coverage",
    tag = "B"
  ) +
  coord_cartesian(clip = "off") +
  theme_paper() +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "plain", color = "black"),
    plot.tag.position = c(-0.05, 1.0)
  )

legend_dat <- tibble::tibble(
  label = unname(lorenz_group_labels[group_order]),
  group = factor(group_order, levels = group_order),
  x = c(0.05, 0.52, 0.05, 0.52),
  y = c(0.68, 0.68, 0.28, 0.28)
)

legend_plot <- ggplot(legend_dat) +
  geom_segment(
    aes(x = x, xend = x + 0.075, y = y, yend = y, color = group),
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_text(
    aes(x = x + 0.095, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 4.0,
    color = "black"
  ) +
  scale_color_manual(values = group_cols, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 8, b = 0, l = 8))

coverage_threshold <- 1.0
bin_file <- "bin_qc_10kb.tsv.gz"

profile_loci <- tibble::tribble(
  ~species, ~locus_label, ~contig, ~feature_start, ~feature_end, ~plot_start, ~plot_end,
  "A. niger", "NT_166526.1:860-870 kb", "NT_166526.1", 860000L, 870000L, 830000L, 900000L,
  "C. nagasakiense", "cntw_is_7:2610-2620 kb", "cntw_is_7", 2610000L, 2620000L, 2580000L, 2650000L
) %>%
  mutate(
    panel = c("i", "ii"),
    locus_key = paste(species, locus_label, sep = "|"),
    title = paste0(panel, "     ", locus_label)
  )

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
        filter(contig == loc$contig, end > loc$plot_start, start < loc$plot_end) %>%
        transmute(
          sample = sample,
          species = species_label,
          method = method_label,
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

an_profile_loci <- profile_loci %>% filter(species == "A. niger")
cntw_profile_loci <- profile_loci %>% filter(species == "C. nagasakiense")

profile_raw <- bind_rows(
  read_bins_for_method(file.path(cluster_root, "AN_shallow"), "Direct route", "A. niger", an_profile_loci),
  read_bins_for_method(file.path(cluster_root, "AN_deep"), "Indirect route", "A. niger", an_profile_loci),
  read_bins_for_method(file.path(cluster_root, "CNTW_shallow"), "Direct route", "C. nagasakiense", cntw_profile_loci),
  read_bins_for_method(file.path(cluster_root, "CNTW_deep"), "Indirect route", "C. nagasakiense", cntw_profile_loci)
)

profile_dat <- profile_raw %>%
  group_by(
    species, method, locus_key, locus_label, panel, title, contig,
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
    colour_key = recode(
      paste(species, method, sep = "|"),
      "A. niger|Direct route" = "AN_direct",
      "A. niger|Indirect route" = "AN_indirect",
      "C. nagasakiense|Direct route" = "CNTW_direct",
      "C. nagasakiense|Indirect route" = "CNTW_indirect"
    )
  )

inner_x_breaks <- function(lims) {
  br <- pretty(lims, n = 4)
  br <- br[br > lims[1] & br < lims[2]]
  if (length(br) < 2) br <- pretty(lims, n = 3)
  br
}

theme_profile <- function(base_size = 8.8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.55),
      axis.ticks = element_line(color = "black", linewidth = 0.55),
      plot.title = element_text(size = base_size + 1.8, face = "plain", color = "black", hjust = 0.02),
      plot.margin = margin(2, 5, 2, 5)
    )
}

make_profile_panel <- function(loc) {
  dat <- profile_dat %>% filter(locus_key == loc$locus_key)
  x_lims <- c(loc$plot_start / 1000, loc$plot_end / 1000)
  dat <- dat %>% filter(mid_kb >= x_lims[1], mid_kb <= x_lims[2])
  highlight_base <- data.frame(xmin = loc$feature_start / 1000, xmax = loc$feature_end / 1000)
  gc_dat <- dat %>% distinct(mid_kb, gc_frac, .keep_all = TRUE) %>% arrange(mid_kb)
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
    theme_profile() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_depth <- ggplot(dat, aes(x = mid_kb, y = mean_norm_depth, color = colour_key, fill = colour_key)) +
    geom_rect(data = highlight_depth, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey90", alpha = 0.35) +
    geom_area(alpha = 0.10, linewidth = 0, position = "identity") +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, color = "grey60", linetype = "dotted", linewidth = 0.35) +
    annotate("text", x = x_lims[2], y = y_max, label = sprintf("[0-%.0fx]", y_max), hjust = 1, vjust = 1, size = 2.3, color = "grey30") +
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    scale_x_continuous(limits = x_lims, breaks = inner_x_breaks(x_lims), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(limits = c(0, depth_lim_hi), breaks = pretty_breaks(n = 3), expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Norm. depth") +
    guides(color = "none", fill = "none") +
    theme_profile() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_cov <- ggplot(dat, aes(x = mid_kb, y = prop_covered, color = colour_key, fill = colour_key)) +
    geom_rect(data = highlight_cov, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey90", alpha = 0.35) +
    geom_area(alpha = 0.10, linewidth = 0, position = "identity") +
    geom_line(linewidth = 0.65) +
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    scale_x_continuous(limits = x_lims, breaks = inner_x_breaks(x_lims), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Position on reference (kb)", y = "Cov.") +
    guides(color = "none", fill = "none") +
    theme_profile()

  p_gc / p_depth / p_cov + plot_layout(heights = c(0.82, 2.35, 0.62))
}

species_title <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4.2, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))
}

profile_panels <- split(profile_loci, seq_len(nrow(profile_loci))) %>% lapply(make_profile_panel)
p_c <- (
  (
    species_title("A. niger") +
      labs(tag = "C") +
      theme(
        plot.tag = element_text(size = 16, face = "plain", color = "black"),
        plot.tag.position = c(0, 0.5)
      ) |
      species_title("C. nagasakiense")
  ) /
    (profile_panels[[1]] | profile_panels[[2]])
) + plot_layout(heights = c(0.12, 1))

top_ab <- (p_gini | p_lorenz) +
  plot_layout(widths = c(1.05, 1))

p_combined <- top_ab / p_c / legend_plot +
  plot_layout(heights = c(1.02, 0.82, 0.14))

ggsave(
  filename = file.path(outdir, "UPDATED_fig4AB_gini_lorenz_memo_publication.png"),
  plot = p_combined,
  width = 12.2,
  height = 9.6,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "UPDATED_fig4AB_gini_lorenz_memo_publication.pdf"),
  plot = p_combined,
  width = 12.2,
  height = 9.6,
  bg = "white",
  limitsize = FALSE
)

write_csv(gini_dat, file.path(outdir, "fig4A_gini_by_group_data.csv"))
write_csv(gini_group_summary, file.path(outdir, "fig4A_gini_by_group_summary.csv"))
write_csv(lorenz_dat, file.path(outdir, "fig4B_lorenz_input_data.csv"))
write_csv(lorenz_summary, file.path(outdir, "fig4B_lorenz_summary_data.csv"))
write_csv(lorenz_group_summary, file.path(outdir, "fig4B_lorenz_group_summary.csv"))

message("Done.")
message("Saved to: ", normalizePath(outdir))
