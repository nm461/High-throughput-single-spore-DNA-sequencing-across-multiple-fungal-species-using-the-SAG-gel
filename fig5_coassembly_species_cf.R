#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
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
if (length(args) >= 1) {
  repo_root <- normalizePath(args[1], mustWork = TRUE)
}

in_dir <- file.path(repo_root, "output", "fig5_coassembly")
out_dir <- file.path(repo_root, "outputs", "UPDATED_fig5_coassembly_species_memo_publication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_root <- dirname(dirname(repo_root))
coassembly_dir <- file.path(data_root, "assembly", "co_assembly")
direct_bootstrap_n <- 40L
direct_bootstrap_paths <- c(
  "A.niger" = file.path(coassembly_dir, "AN_shallow_direct_bootstrap_BUSCO_summary.csv"),
  "C.nagasakiense" = file.path(coassembly_dir, "CNTW_shallow_direct_bootstrap_BUSCO_summary.csv")
)
indirect_path <- file.path(in_dir, "fig5_MDAx2_sequential_input.csv")

indirect_raw <- readr::read_csv(indirect_path, show_col_types = FALSE)

species_levels <- c("A.niger", "C.nagasakiense")
species_labels <- c(
  "A.niger" = "A. niger",
  "C.nagasakiense" = "C. nagasakiense"
)

metric_levels <- c("Completeness", "Fragmentation", "Missing")

route_cols <- c(
  "A.niger|Direct route" = "#E1BE6A",
  "A.niger|Indirect route" = "#1A85FF",
  "C.nagasakiense|Direct route" = "#2D8B82",
  "C.nagasakiense|Indirect route" = "#D41159"
)

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x = element_text(size = base_size + 1, color = "black", margin = margin(t = 8)),
      axis.title.y = element_text(size = base_size + 1, color = "black", margin = margin(r = 8)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1, face = "plain", color = "black", margin = margin(2, 0, 2, 0)),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.13, "cm"),
      plot.margin = margin(8, 12, 10, 12),
      plot.tag = element_text(size = base_size + 4, face = "plain", color = "black")
    )
}

format_species <- function(x) {
  factor(x, levels = species_levels, labels = unname(species_labels[species_levels]))
}

summarise_direct_bootstrap <- function(path, species_id, n_bootstrap = direct_bootstrap_n) {
  if (!file.exists(path)) stop("Missing direct-route bootstrap input: ", path)

  raw <- readr::read_csv(path, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    mutate(
      iteration = as.integer(iteration),
      sag = as.integer(sag)
    ) %>%
    filter(!is.na(iteration))

  keep_iterations <- raw %>%
    distinct(iteration) %>%
    arrange(iteration) %>%
    slice_head(n = n_bootstrap) %>%
    pull(iteration)

  raw <- raw %>%
    filter(iteration %in% keep_iterations)

  observed_n <- raw %>%
    distinct(iteration) %>%
    nrow()

  if (observed_n != n_bootstrap) {
    stop(sprintf(
      "%s has %d bootstrap iterations after filtering, expected %d",
      basename(path), observed_n, n_bootstrap
    ))
  }

  raw %>%
    select(iteration, sag, completeness, fragmentation, missing) %>%
    tidyr::pivot_longer(
      cols = c(completeness, fragmentation, missing),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      species = species_id,
      route = "MDAx1",
      metric = recode(
        metric,
        completeness = "Completeness",
        fragmentation = "Fragmentation",
        missing = "Missing"
      )
    ) %>%
    group_by(species, route, sag, metric) %>%
    summarise(
      n = n_distinct(iteration),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = paste(route, species))
}

direct_raw <- bind_rows(lapply(names(direct_bootstrap_paths), function(species_id) {
  summarise_direct_bootstrap(direct_bootstrap_paths[[species_id]], species_id)
}))

direct_count_check <- direct_raw %>%
  distinct(species, sag, metric, n) %>%
  count(species, n, name = "rows")

if (!all(direct_count_check$n == direct_bootstrap_n)) {
  stop("Direct-route bootstrap summaries do not all use n=", direct_bootstrap_n)
}

direct <- direct_raw %>%
  transmute(
    species,
    route = "Direct route",
    sag = as.integer(sag),
    metric = factor(metric, levels = metric_levels),
    value = as.numeric(median),
    q25 = as.numeric(q25),
    q75 = as.numeric(q75),
    n_bootstrap = as.integer(n),
    source = paste0("Bootstrap median (first ", direct_bootstrap_n, " iterations)")
  )

indirect <- indirect_raw %>%
  transmute(
    species,
    route = "Indirect route",
    sag = as.integer(sag),
    metric = factor(metric, levels = metric_levels),
    value = as.numeric(value),
    q25 = NA_real_,
    q75 = NA_real_,
    n_bootstrap = NA_integer_,
    source = "Top-ranked sequential"
  )

plot_dat <- bind_rows(direct, indirect) %>%
  filter(species %in% species_levels) %>%
  mutate(
    species_label = format_species(species),
    route = factor(route, levels = c("Direct route", "Indirect route")),
    colour_key = paste(species, route, sep = "|")
  )

completion <- plot_dat %>%
  filter(metric == "Completeness")

endpoint_labels <- completion %>%
  group_by(species, route) %>%
  filter(sag == max(sag, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    label = sprintf("%.1f%%", value),
    label_x = sag + 0.16,
    label_y = case_when(
      route == "Direct route" & species == "A.niger" ~ value + 0.5,
      route == "Direct route" & species == "C.nagasakiense" ~ value - 2.4,
      route == "Indirect route" & species == "A.niger" ~ value - 3.0,
      TRUE ~ value + 1.4
    ),
    hjust_lab = 0
  )

panel_note <- completion %>%
  filter(route == "Direct route") %>%
  group_by(species) %>%
  summarise(n_bootstrap = max(n_bootstrap, na.rm = TRUE), .groups = "drop") %>%
  mutate(note = paste0(n_bootstrap, " bootstrap co-assemblies per SAG count"))

make_completion_panel <- function(sp, tag) {
  dat_sp <- completion %>% filter(species == sp)
  endpoints_sp <- endpoint_labels %>% filter(species == sp)
  note_sp <- panel_note %>% filter(species == sp)
  cols <- route_cols[paste(sp, c("Direct route", "Indirect route"), sep = "|")]
  names(cols) <- c("Direct route", "Indirect route")

  ggplot(dat_sp, aes(x = sag, y = value, color = route, fill = route)) +
    geom_ribbon(
      data = dat_sp %>% filter(route == "Direct route"),
      aes(ymin = q25, ymax = q75),
      alpha = 0.28,
      linewidth = 0,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 2.25) +
    geom_text(
      data = endpoints_sp,
      aes(x = label_x, y = label_y, label = label, color = route),
      inherit.aes = FALSE,
      hjust = endpoints_sp$hjust_lab,
      size = 3.85,
      fontface = "plain"
    ) +
    facet_wrap(~ species_label, nrow = 1) +
    scale_color_manual(values = cols, breaks = c("Direct route", "Indirect route")) +
    scale_fill_manual(values = cols, breaks = c("Direct route", "Indirect route")) +
    scale_x_continuous(breaks = 2:12, limits = c(1.7, 13.35), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_percent(scale = 1), breaks = seq(0, 100, 20), limits = c(0, 105), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Number of SAGs co-assembled", y = "BUSCO completeness") +
    guides(color = "none", fill = "none") +
    theme_pub(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.tag.position = c(0.01, 0.99)
    ) +
    labs(tag = tag)
}

p_an <- make_completion_panel("A.niger", "A")
p_cn <- make_completion_panel("C.nagasakiense", "B")

species_legend_plot <- function() {
  legend_dat <- tibble::tibble(
    label = c(
      "Direct route: A. niger",
      "Indirect route: A. niger",
      "Direct route: C. nagasakiense",
      "Indirect route: C. nagasakiense"
    ),
    colour_key = c(
      "A.niger|Direct route",
      "A.niger|Indirect route",
      "C.nagasakiense|Direct route",
      "C.nagasakiense|Indirect route"
    ),
    x = c(0.06, 0.06, 0.52, 0.52),
    y = c(0.70, 0.23, 0.70, 0.23),
    text_x = c(0.155, 0.155, 0.615, 0.615)
  )

  ggplot(legend_dat) +
    geom_segment(aes(x = x, xend = x + 0.075, y = y, yend = y, color = colour_key), linewidth = 1.2, lineend = "round") +
    geom_point(aes(x = x + 0.0375, y = y, color = colour_key), size = 2.4) +
    geom_text(aes(x = text_x, y = y, label = label), hjust = 0, vjust = 0.5, size = 4.1, color = "black") +
    scale_color_manual(values = route_cols, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 4, b = 0, l = 4))
}

main_fig <- (p_an + p_cn) / species_legend_plot() +
  plot_layout(ncol = 1, heights = c(1, 0.15))

ggsave(
  file.path(out_dir, "UPDATED_fig5_coassembly_species_memo_publication.png"),
  main_fig,
  width = 12.0,
  height = 6.45,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "UPDATED_fig5_coassembly_species_memo_publication.pdf"),
  main_fig,
  width = 12.0,
  height = 6.45,
  bg = "white",
  limitsize = FALSE
)

write_csv(plot_dat, file.path(out_dir, "UPDATED_fig5_coassembly_species_memo_plot_data.csv"))
write_csv(endpoint_labels, file.path(out_dir, "UPDATED_fig5_coassembly_species_memo_endpoint_labels.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir))

# ── By-species figure with Completeness + Fragmentation strip titles ──────────

cf_dat <- plot_dat %>%
  filter(metric %in% c("Completeness", "Fragmentation")) %>%
  mutate(metric = factor(metric, levels = c("Completeness", "Fragmentation")))

cf_endpoint_labels <- cf_dat %>%
  group_by(route, species, metric) %>%
  filter(sag == max(sag, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    label = sprintf("%.1f%%", value),
    label_x = sag + 0.16,
    label_y = case_when(
      metric == "Completeness" & route == "Direct route"   & species == "A.niger"         ~ value + 0.5,
      metric == "Completeness" & route == "Direct route"   & species == "C.nagasakiense"  ~ value - 2.4,
      metric == "Completeness" & route == "Indirect route" & species == "A.niger"         ~ value - 3.0,
      metric == "Completeness"                                                              ~ value + 1.4,
      metric == "Fragmentation" & route == "Direct route"  & species == "A.niger"         ~ value + 0.8,
      metric == "Fragmentation" & route == "Direct route"  & species == "C.nagasakiense"  ~ value - 0.9,
      metric == "Fragmentation" & route == "Indirect route"& species == "A.niger"         ~ value - 1.0,
      TRUE                                                                                  ~ value - 0.9
    ),
    hjust_lab = 0
  )

metric_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 5.0, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_cf_species_panel <- function(sp, metric_name, tag,
                                   show_x_axis = FALSE, show_y_axis = TRUE,
                                   show_y_ticks = show_y_axis) {
  dat    <- cf_dat %>% filter(species == sp, metric == metric_name)
  labels <- cf_endpoint_labels %>% filter(species == sp, metric == metric_name)
  cols   <- route_cols[paste(sp, c("Direct route", "Indirect route"), sep = "|")]
  names(cols) <- c("Direct route", "Indirect route")

  y_limits <- if (metric_name == "Completeness") c(0, 100) else c(0, 20)
  y_breaks <- if (metric_name == "Completeness") seq(0, 100, 25) else seq(0, 20, 5)
  y_title  <- if (show_y_axis) paste0("BUSCO ", tolower(metric_name)) else NULL
  sp_label <- species_labels[[sp]]

  p <- ggplot(dat, aes(x = sag, y = value, color = route, fill = route)) +
    geom_ribbon(
      data = dat %>% filter(route == "Direct route"),
      aes(ymin = q25, ymax = q75),
      alpha = 0.28, linewidth = 0, show.legend = FALSE
    ) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.0) +
    geom_text(
      data = labels,
      aes(x = label_x, y = label_y, label = label, color = route),
      inherit.aes = FALSE,
      hjust = 0, size = 3.35, fontface = "plain"
    ) +
    facet_wrap(~ factor(sp_label, levels = sp_label), nrow = 1) +
    scale_color_manual(values = cols, breaks = c("Direct route", "Indirect route")) +
    scale_fill_manual(values  = cols, breaks = c("Direct route", "Indirect route")) +
    scale_x_continuous(breaks = 2:12, limits = c(1.7, 13.6),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_percent(scale = 1),
                       breaks = y_breaks, limits = y_limits,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      x   = if (show_x_axis) "Number of SAGs co-assembled" else NULL,
      y   = y_title,
      tag = tag
    ) +
    guides(color = "none", fill = "none") +
    theme_pub(base_size = 12) +
    theme(
      axis.title.x  = if (show_x_axis) element_text(size = 12.5, margin = margin(t = 8)) else element_blank(),
      axis.text.x   = if (show_x_axis) element_text(size = 10.5) else element_blank(),
      axis.ticks.x  = if (show_x_axis) element_line(color = "black", linewidth = 0.8) else element_blank(),
      axis.title.y  = if (show_y_axis) element_text(size = 12.5, margin = margin(r = 8)) else element_blank(),
      axis.text.y   = if (show_y_ticks) element_text(size = 10.5) else element_blank(),
      axis.ticks.y  = if (show_y_ticks) element_line(color = "black", linewidth = 0.8) else element_blank(),
      plot.tag      = element_text(size = 13, face = "plain"),
      plot.tag.position = c(0.01, 0.99),
      plot.margin   = margin(4, 8, 4, 8)
    )
  p
}

# Row panels
comp_an   <- make_cf_species_panel("A.niger",        "Completeness",  "A", show_x_axis = TRUE, show_y_axis = TRUE)
comp_cn   <- make_cf_species_panel("C.nagasakiense", "Completeness",  "",  show_x_axis = TRUE, show_y_axis = FALSE, show_y_ticks = TRUE)
frag_an   <- make_cf_species_panel("A.niger",        "Fragmentation", "B", show_x_axis = TRUE,  show_y_axis = TRUE)
frag_cn   <- make_cf_species_panel("C.nagasakiense", "Fragmentation", "",  show_x_axis = TRUE,  show_y_axis = FALSE, show_y_ticks = TRUE)

completeness_row   <- comp_an + comp_cn + plot_layout(ncol = 2)
fragmentation_row  <- frag_an + frag_cn + plot_layout(ncol = 2)

species_cf_fig <- metric_strip("Completeness") / completeness_row /
  patchwork::plot_spacer() /
  metric_strip("Fragmentation") / fragmentation_row /
  species_legend_plot() +
  plot_layout(ncol = 1, heights = c(0.125, 1, 0.06, 0.125, 1, 0.15))

ggsave(
  file.path(out_dir, "fig5_coassembly_species_completeness_fragmentation.png"),
  species_cf_fig,
  width = 12.0, height = 9.6, dpi = 300, bg = "white", limitsize = FALSE
)
ggsave(
  file.path(out_dir, "fig5_coassembly_species_completeness_fragmentation.pdf"),
  species_cf_fig,
  width = 12.0, height = 9.6, bg = "white", limitsize = FALSE
)

message("Species CF figure saved.")

metrics_dat <- plot_dat %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    route_species = factor(
      paste(route, species_label, sep = ": "),
      levels = c(
        "Direct route: A. niger",
        "Indirect route: A. niger",
        "Direct route: C. nagasakiense",
        "Indirect route: C. nagasakiense"
      )
    )
  )

diagnostic_cols <- c(
  "Direct route: A. niger" = route_cols[["A.niger|Direct route"]],
  "Indirect route: A. niger" = route_cols[["A.niger|Indirect route"]],
  "Direct route: C. nagasakiense" = route_cols[["C.nagasakiense|Direct route"]],
  "Indirect route: C. nagasakiense" = route_cols[["C.nagasakiense|Indirect route"]]
)

diagnostic_fig <- ggplot(metrics_dat, aes(x = sag, y = value, color = route_species, fill = route_species)) +
  geom_ribbon(
    data = metrics_dat %>% filter(route == "Direct route"),
    aes(ymin = q25, ymax = q75),
    alpha = 0.24,
    linewidth = 0,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.85) +
  facet_grid(metric ~ species_label, scales = "free_y") +
  scale_color_manual(values = diagnostic_cols) +
  scale_fill_manual(values = diagnostic_cols) +
  scale_x_continuous(breaks = 2:12, limits = c(1.7, 12.3)) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "Number of SAGs co-assembled", y = "BUSCO metric") +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.3, size = 2.8, alpha = 1)),
    fill = "none"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(out_dir, "fig5_coassembly_busco_metrics_publication.png"),
  diagnostic_fig,
  width = 10.8,
  height = 8.6,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "fig5_coassembly_busco_metrics_publication.pdf"),
  diagnostic_fig,
  width = 10.8,
  height = 8.6,
  bg = "white",
  limitsize = FALSE
)

write_csv(plot_dat, file.path(out_dir, "fig5_coassembly_publication_plot_data.csv"))
write_csv(endpoint_labels, file.path(out_dir, "fig5_coassembly_publication_endpoint_labels.csv"))

route_endpoint_labels <- completion %>%
  group_by(route, species) %>%
  filter(sag == max(sag, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    label = sprintf("%.1f%%", value),
    label_x = sag + 0.16,
    label_y = case_when(
      route == "Direct route" & species == "A.niger" ~ value + 0.5,
      route == "Direct route" & species == "C.nagasakiense" ~ value - 2.4,
      route == "Indirect route" & species == "A.niger" ~ value - 3.0,
      TRUE ~ value + 1.4
    ),
    hjust_lab = 0
  )

make_route_legend <- function() {
  legend_dat <- tibble(
    label = c(
      "Direct route: A. niger",
      "Direct route: C. nagasakiense",
      "Indirect route: A. niger",
      "Indirect route: C. nagasakiense"
    ),
    colour_key = c(
      "A.niger|Direct route",
      "C.nagasakiense|Direct route",
      "A.niger|Indirect route",
      "C.nagasakiense|Indirect route"
    ),
    x = c(0.15, 0.56, 0.15, 0.56),
    y = c(0.68, 0.68, 0.28, 0.28)
  )

  ggplot(legend_dat) +
    geom_segment(
      aes(x = x, xend = x + 0.075, y = y, yend = y, color = colour_key),
      linewidth = 1.25,
      lineend = "round"
    ) +
    geom_point(aes(x = x + 0.0375, y = y, color = colour_key), size = 2.4) +
    geom_text(aes(x = x + 0.095, y = y, label = label), hjust = 0, vjust = 0.5, size = 4.2, color = "black") +
    scale_color_manual(values = route_cols, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 4, b = 0, l = 4))
}

method_note_plot <- function() {
  route_note <- paste0(
    "Direct route: bootstrap median with shaded IQR; A. niger n=53 and C. nagasakiense n=44 bootstrap co-assemblies per SAG count.\n",
    "Indirect route: top-ranked sequential co-assembly."
  )

  ggplot() +
    geom_text(aes(x = 0.02, y = 0.58, label = route_note), size = 3.35, color = "black", hjust = 0, lineheight = 1.02) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
}

make_route_panel <- function(route_name, tag) {
  dat_route <- completion %>% filter(route == route_name)
  endpoints_route <- route_endpoint_labels %>% filter(route == route_name)
  cols <- route_cols[paste(species_levels, route_name, sep = "|")]

  ggplot(dat_route, aes(x = sag, y = value, color = colour_key, fill = colour_key)) +
    geom_ribbon(
      data = dat_route %>% filter(route == "Direct route"),
      aes(ymin = q25, ymax = q75),
      alpha = 0.34,
      linewidth = 0,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 2.25) +
    geom_text(
      data = endpoints_route,
      aes(x = label_x, y = label_y, label = label, color = colour_key),
      inherit.aes = FALSE,
      hjust = endpoints_route$hjust_lab,
      size = 3.9,
      fontface = "plain"
    ) +
    facet_wrap(~ route, nrow = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = 2:12, limits = c(1.7, 13.6), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_percent(scale = 1), breaks = seq(0, 100, 20), limits = c(0, 105), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Number of SAGs co-assembled", y = "BUSCO completeness") +
    guides(color = "none", fill = "none") +
    theme_pub(base_size = 13) +
    theme(plot.tag.position = c(0.01, 0.99)) +
    labs(tag = tag)
}

p_direct_route <- make_route_panel("Direct route", "A")
p_indirect_route <- make_route_panel("Indirect route", "B")

route_grouped_fig <- (p_direct_route + p_indirect_route) / make_route_legend() / method_note_plot() +
  plot_layout(ncol = 1, heights = c(1, 0.13, 0.15))

ggsave(
  file.path(out_dir, "fig5_coassembly_completeness_by_route_publication.png"),
  route_grouped_fig,
  width = 12.0,
  height = 6.15,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "fig5_coassembly_completeness_by_route_publication.pdf"),
  route_grouped_fig,
  width = 12.0,
  height = 6.15,
  bg = "white",
  limitsize = FALSE
)

cf_dat <- plot_dat %>%
  filter(metric %in% c("Completeness", "Fragmentation")) %>%
  mutate(
    metric = factor(metric, levels = c("Completeness", "Fragmentation")),
    colour_key = paste(species, route, sep = "|")
  )

cf_endpoint_labels <- cf_dat %>%
  group_by(route, species, metric) %>%
  filter(sag == max(sag, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    label = sprintf("%.1f%%", value),
    label_x = sag + 0.16,
    label_y = case_when(
      metric == "Completeness" & route == "Direct route" & species == "A.niger" ~ value + 0.5,
      metric == "Completeness" & route == "Direct route" & species == "C.nagasakiense" ~ value - 2.4,
      metric == "Completeness" & route == "Indirect route" & species == "A.niger" ~ value - 3.0,
      metric == "Completeness" ~ value + 1.4,
      metric == "Fragmentation" & route == "Direct route" & species == "A.niger" ~ value + 0.8,
      metric == "Fragmentation" & route == "Direct route" & species == "C.nagasakiense" ~ value - 0.9,
      metric == "Fragmentation" & route == "Indirect route" & species == "A.niger" ~ value - 1.0,
      TRUE ~ value - 0.9
    ),
    hjust_lab = ifelse(label_x > sag, 0, 1)
  )

metric_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 5.0, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_cf_panel <- function(route_name, metric_name, tag, show_x_axis = FALSE, show_y_axis = TRUE, show_y_ticks = show_y_axis) {
  dat <- cf_dat %>% filter(route == route_name, metric == metric_name)
  labels <- cf_endpoint_labels %>% filter(route == route_name, metric == metric_name)
  cols <- route_cols[paste(species_levels, route_name, sep = "|")]

  y_limits <- if (metric_name == "Completeness") c(0, 100) else c(0, 20)
  y_breaks <- if (metric_name == "Completeness") seq(0, 100, 25) else seq(0, 20, 5)
  y_title <- if (show_y_axis) paste0("BUSCO ", tolower(metric_name)) else NULL

  ggplot(dat, aes(x = sag, y = value, color = colour_key, fill = colour_key)) +
    geom_ribbon(
      data = dat %>% filter(route == "Direct route"),
      aes(ymin = q25, ymax = q75),
      alpha = 0.34,
      linewidth = 0,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.0) +
    geom_text(
      data = labels,
      aes(x = label_x, y = label_y, label = label, color = colour_key),
      inherit.aes = FALSE,
      hjust = labels$hjust_lab,
      size = 3.35,
      fontface = "plain"
    ) +
    facet_wrap(~ route, nrow = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = 2:12, limits = c(1.7, 13.6), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_percent(scale = 1), breaks = y_breaks, limits = y_limits, expand = expansion(mult = c(0, 0.05))) +
    labs(x = if (show_x_axis) "Number of SAGs co-assembled" else NULL, y = y_title, tag = tag) +
    guides(color = "none", fill = "none") +
    theme_pub(base_size = 12) +
    theme(
      axis.title.x = if (show_x_axis) element_text(size = 12.5, margin = margin(t = 8)) else element_blank(),
      axis.text.x = if (show_x_axis) element_text(size = 10.5) else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line(color = "black", linewidth = 0.8) else element_blank(),
      axis.title.y = if (show_y_axis) element_text(size = 12.5, margin = margin(r = 8)) else element_blank(),
      axis.text.y = if (show_y_ticks) element_text(size = 10.5) else element_blank(),
      axis.ticks.y = if (show_y_ticks) element_line(color = "black", linewidth = 0.8) else element_blank(),
      plot.tag = element_text(size = 13, face = "plain"),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(4, 8, 4, 8)
    )
}

cf_completeness_row <- make_cf_panel("Direct route", "Completeness", "A", show_x_axis = FALSE, show_y_axis = TRUE) +
  make_cf_panel("Indirect route", "Completeness", "", show_x_axis = FALSE, show_y_axis = FALSE, show_y_ticks = TRUE) +
  plot_layout(ncol = 2)

cf_fragmentation_row <- make_cf_panel("Direct route", "Fragmentation", "B", show_x_axis = TRUE, show_y_axis = TRUE) +
  make_cf_panel("Indirect route", "Fragmentation", "", show_x_axis = TRUE, show_y_axis = FALSE, show_y_ticks = TRUE) +
  plot_layout(ncol = 2)

cf_route_fig <- metric_strip("Completeness") / cf_completeness_row /
  patchwork::plot_spacer() /
  metric_strip("Fragmentation") / cf_fragmentation_row /
  make_route_legend() / method_note_plot() +
  plot_layout(ncol = 1, heights = c(0.125, 1, 0.08, 0.125, 1, 0.15, 0.15))

ggsave(
  file.path(out_dir, "fig5_coassembly_completeness_fragmentation_by_route_publication.png"),
  cf_route_fig,
  width = 11.0,
  height = 9.4,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "fig5_coassembly_completeness_fragmentation_by_route_publication.pdf"),
  cf_route_fig,
  width = 11.0,
  height = 9.4,
  bg = "white",
  limitsize = FALSE
)

message("Done. Outputs written to: ", normalizePath(out_dir))
