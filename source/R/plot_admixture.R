library(tidyverse)
library(patchwork)
library(ggdendro)

args <- commandArgs(trailing = TRUE)

k_min_error <- as.integer(args[2])

sample_metadata <- read.table(
  args[3],
  sep = ",",
  header = FALSE,
  col.names = c("ID", "Species", "Population", "Sex")
)

# Only the first two columns are named as additional metadata may be added later
read_group_metadata <- function(path) {
  group_metadata <- read.table(
    path,
    sep = ",",
    header = FALSE,
    comment.char = "" # Avoids hex codes reading in as comments
  )
  names(group_metadata)[1:2] <- c("GROUP", "COLOUR")
  group_metadata
}

group_palette <- function(group_metadata) {
  setNames(group_metadata$COLOUR, group_metadata$GROUP)
}

population_palette <- group_palette(read_group_metadata(args[4]))
species_palette <- group_palette(read_group_metadata(args[5]))

admixture <- data.table::fread(args[1], fill = Inf) |>
  as_tibble() |>
  rename(K = V1, ID = V2)

admixture_long <- admixture |>
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "CLUSTER",
    names_transform = \(x) as.integer(x) - 2,
    values_to = "PROBABILITY",
    values_drop_na = TRUE
  )

admixture_clusters <- admixture |>
  filter(K == max(K)) |>
  column_to_rownames("ID") |>
  select(!K) |>
  dist() |>
  hclust(method = "ward.D2")

dendro_data <- as.dendrogram(admixture_clusters) |> dendro_data()
ids_sorted <- data.frame(
  ID = dendro_data$labels$label,
  SEQUENCE = dendro_data$labels$x
)

admixture_long <- admixture_long |> left_join(ids_sorted)
xmin <- min(admixture_long$SEQUENCE)
xmax <- max(admixture_long$SEQUENCE)

sample_metadata_sorted <- admixture_long |>
  left_join(sample_metadata) |>
  arrange(SEQUENCE) |>
  mutate(
    Population = factor(Population, levels = unique(Population)),
    Species = factor(Species, levels = unique(Species))
  )

admixture_plot <- admixture_long |>
  ggplot(
    aes(
      x = SEQUENCE,
      y = PROBABILITY,
      fill = as.factor(CLUSTER),
      alpha = (K == k_min_error)
    )
  ) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  geom_col(width = 1, show.legend = FALSE) +
  facet_grid(
    rows = vars(K),
    switch = "y",
    labeller = label_bquote(K == .(K))
  ) +
  scale_alpha_manual(values = c(0.33, 1)) +
  scale_fill_brewer(palette = "Set1") +
  theme_void() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.15),
    strip.clip = "off",
    strip.text.y.left = element_text(
      hjust = 0,
      vjust = 1,
      angle = 0,
      size = 6,
      margin = margin(t = 1.5, r = -7.5, unit = "mm")
    )
  )

dendrogram <- ggdendro::segment(dendro_data) |>
  ggplot(aes(x = x, y = -y, xend = xend, yend = -yend)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  geom_segment(linewidth = 0.15) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm")
  )

pop_meta <- sample_metadata_sorted |>
  ggplot(aes(x = SEQUENCE, y = 0, fill = Population)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  geom_tile() +
  scale_fill_manual(values = population_palette, drop = TRUE) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(colour = "black", linewidth = 0.15)
    )
  ) +
  theme_void() +
  theme(
    plot.margin = unit(c(2, 0, 0, 0), units = "mm"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  )

spp_meta <- sample_metadata_sorted |>
  ggplot(aes(x = SEQUENCE, y = 0, fill = Species)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  geom_tile() +
  scale_fill_manual(values = species_palette, drop = TRUE) +
  guides(
    fill = guide_legend(
      order = 2,
      override.aes = list(colour = "black", linewidth = 0.15)
    )
  ) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  )

combined_plot <- (admixture_plot / pop_meta / dendrogram) +
  plot_layout(
    guides = "collect",
    heights = c(
      43.5 - (length(unique(admixture$K))) / 2,
      0 + (length(unique(admixture$K))) / 2,
      6.5
    )
  )

if (length(unique(sample_metadata$Population)) > 1 && length(unique(sample_metadata$Species)) > 1) {
  combined_plot <- (admixture_plot / pop_meta / spp_meta / dendrogram) +
    plot_layout(
      guides = "collect",
      heights = c(
        44 - (length(unique(admixture$K))) / 2,
        0 + (length(unique(admixture$K))) / 2,
        0 + (length(unique(admixture$K))) / 2,
        6
      )
    )
}

combined_plot <- combined_plot &
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.key.size = unit(2, "mm"),
    legend.key.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(t = 0, b = -1, unit = "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold")
  )

small_plot <- combined_plot &
  theme(
    legend.position = "none",
    strip.text.y.left = element_text(
      hjust = 0,
      vjust = 1,
      angle = 0,
      size = 7,
      margin = margin(t = 1.5, r = -7.5, unit = "mm")
    )
  )

ggsave(
  plot = combined_plot,
  file = "admixture.png",
  dpi = 600,
  width = 6.75,
  height = 1 + (0.3 * length(unique(admixture$K))),
  bg = "white"
)

ggsave(
  plot = small_plot,
  file = "admixture_small.png",
  dpi = 600,
  width = (6.75 / 2),
  height = 1 + (0.15 * length(unique(admixture$K))),
  bg = "white"
)
