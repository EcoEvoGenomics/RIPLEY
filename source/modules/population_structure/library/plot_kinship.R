library(tidyverse)
library(patchwork)
library(ggdendro)
library(scico)

args <- commandArgs(trailing = TRUE)
name <- basename(args[1])

sample_metadata <- read.table(
  args[2],
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

# Metadata annotation is only informative if the data has > 1 levels
is_informative <- function(x) n_distinct(x, na.rm = TRUE) > 1

population_palette <- group_palette(read_group_metadata(args[3]))
species_palette <- group_palette(read_group_metadata(args[4]))

data <- read.table(args[1], header = TRUE) |>
  select(INDV1, INDV2, RELATEDNESS_PHI) |>
  rename(ID = INDV1, ID2 = INDV2, PHI = RELATEDNESS_PHI) |>
  left_join(sample_metadata)

dist <- data |>
  select(ID, ID2, PHI) |>
  pivot_wider(names_from = ID2, values_from = PHI) |>
  column_to_rownames("ID") |>
  as.matrix()

dist <- dist[rownames(dist), rownames(dist)]
clust <- hclust(as.dist(1 - dist), method = "ward.D2")
dendro_data <- as.dendrogram(clust) |> dendro_data()
ids_sorted <- data.frame(
  ID = dendro_data$labels$label,
  SEQUENCE = dendro_data$labels$x
)

data <- data |>
  left_join(ids_sorted) |>
  left_join(
    ids_sorted |> rename(ID2 = ID, SEQUENCE2 = SEQUENCE),
    by = "ID2"
  ) |>
  filter(SEQUENCE >= SEQUENCE2) |>
  mutate(PHI = ifelse(ID == ID2, NA, PHI))

clustered_levels <- function(data, column) {
  data |>
    arrange(SEQUENCE) |>
    pull({{ column }}) |>
    unique()
}

data <- data |>
  mutate(
    Population = factor(Population, levels = clustered_levels(data, Population)),
    Species = factor(Species, levels = clustered_levels(data, Species))
  )

xmin <- min(data$SEQUENCE) - 0.5
xmax <- max(data$SEQUENCE) + 0.5
phi_bound <- max(abs(data$PHI), na.rm = TRUE)
phi_limits <- c(-phi_bound, phi_bound)

kinship_matrix <- data |>
  ggplot(
    aes(
      x = SEQUENCE,
      y = factor(ID2, levels = ids_sorted$ID),
      fill = PHI
    )
  ) +
  coord_equal(expand = FALSE, xlim = c(xmin, xmax)) +
  scale_fill_scico(
    palette = "vik",
    midpoint = 0,
    limits = phi_limits,
    na.value = "white"
  ) +
  geom_tile(show.legend = FALSE) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 1, 1.5, 1), unit = "mm")
  )

kinship_key <- data.frame(
  PHI = seq(phi_limits[1], phi_limits[2], length.out = nrow(data))
) |>
  ggplot(aes(x = 0, y = PHI, fill = PHI)) +
  ggtitle("Kinship") +
  coord_cartesian(expand = FALSE) +
  geom_raster(show.legend = FALSE) +
  scale_y_continuous(
    position = "left",
    labels = scales::number_format(accuracy = 0.001)
  ) +
  scale_fill_scico(
    palette = "vik",
    midpoint = 0,
    limits = phi_limits
  ) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 5, margin = margin(r = 1, unit = "mm")),
    axis.ticks.length = unit(0.25, units = "mm"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.15),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.15),
    # Bottom margin trims the bar's lower edge flush with the matrix bottom
    plot.margin = unit(c(0, 0, 0.65, 0), unit = "mm"),
    plot.title = element_text(size = 6, face = "bold", hjust = 1, vjust = 4)
  )

dendrogram <- ggdendro::segment(dendro_data) |>
  ggplot(aes(x = x, y = -y, xend = xend, yend = -yend)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  geom_segment(linewidth = 0.15) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm")
  )

pop_meta <- data |>
  ggplot(aes(x = SEQUENCE, y = 0, fill = Population)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  scale_fill_manual(
    values = population_palette,
    drop = TRUE,
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.10)
    )
  ) +
  geom_tile() +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  )

spp_meta <- data |>
  ggplot(aes(x = SEQUENCE, y = 0, fill = Species)) +
  coord_cartesian(expand = FALSE, xlim = c(xmin, xmax)) +
  scale_fill_manual(
    values = species_palette,
    drop = TRUE,
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.10)
    )
  ) +
  geom_tile() +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  )

meta_rows <- c(
  if (is_informative(data$Population)) list(pop_meta),
  if (is_informative(data$Species)) list(spp_meta)
)

matrix_column <- wrap_plots(
  c(list(kinship_matrix), meta_rows, list(dendrogram)),
  ncol = 1,
  guides = "collect",
  heights = c(37, rep(0.5, length(meta_rows)), 3 - (length(meta_rows) / 2))
)

combined_plot <- (
  (plot_spacer() / kinship_key / plot_spacer()) +
    plot_layout(heights = c(0.01, 36.5, 3.5)) |
    matrix_column
) + plot_layout(widths = c(1, 39))

combined_plot <- combined_plot &
  theme(
    legend.position = "left",
    legend.justification = "top",
    legend.box.spacing = unit(0, "mm"),
    legend.key.size = unit(2, "mm"),
    legend.key.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(t = 2.27, b = -2.5, l = -2.4, r = -20.6, unit = "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(
      size = 6,
      face = "bold",
      margin = margin(b = 1.55, unit = "mm")
    )
  )

ggsave(
  plot = combined_plot,
  filename = paste0(name, ".png"),
  dpi = 600,
  width  = 6.75 / 2,
  height = 6.75 / 2,
  bg = "white"
)
