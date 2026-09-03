library(tidyverse)
library(patchwork)
library(ggdendro)

args <- commandArgs(trailing = TRUE)
name <- basename(args[1])

meta <- read.table(
  args[2],
  sep = ",",
  header = FALSE,
  col.names = c("ID", "Species", "Population", "Sex")
)

data <- read.table(args[1], header = TRUE) |>
  select(INDV1, INDV2, RELATEDNESS_PHI) |>
  rename(ID = INDV1, ID2 = INDV2, PHI = RELATEDNESS_PHI) |>
  left_join(meta)

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

xmin <- min(data$SEQUENCE) - 0.5
xmax <- max(data$SEQUENCE) + 0.5

kinship_matrix <- data |>
  ggplot(
    aes(
      x = SEQUENCE,
      y = factor(ID2, levels = ids_sorted$ID),
      fill = PHI
    )
  ) +
  coord_equal(expand = FALSE, xlim = c(xmin, xmax)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white"
  ) +
  geom_tile(show.legend = FALSE) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 1, 1.5, 1), unit = "mm")
  )

kinship_key <- data |>
  ggplot(
    aes(
      x = 0,
      y = seq(min(PHI, na.rm = TRUE), max(PHI, na.rm = TRUE), length.out = length(PHI)),
      fill = seq(min(PHI, na.rm = TRUE), max(PHI, na.rm = TRUE), length.out = length(PHI))
    )
  ) +
  ggtitle(expression(bold("Kinship" ~ Phi))) +
  coord_cartesian(expand = FALSE) +
  geom_raster(show.legend = FALSE) +
  scale_y_continuous(
    position = "right",
    labels = scales::number_format(accuracy = 0.001)
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 5, margin = margin(l = 1, unit = "mm")),
    axis.ticks.length = unit(0.25, units = "mm"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.15),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.15),
    plot.margin = unit(c(0, 0, 1, 0), unit = "mm"),
    plot.title = element_text(size = 6, face = "bold", hjust = 0, vjust = 4)
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
  scale_fill_discrete(
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.15)
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
  scale_fill_discrete(
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.15)
    )
  ) +
  geom_tile() +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "mm"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  )

combined_plot <- (
  (kinship_matrix / pop_meta / dendrogram) +
    plot_layout(guides = "collect", heights = c(37, 0.5, 2.5)
    ) |
    (plot_spacer() / kinship_key / plot_spacer()) +
      plot_layout(heights = c(0.01, 36.6, 3.4)
      )
) + plot_layout(widths = c(39, 1))

if (length(unique(data$Population)) > 1 && length(unique(data$Species)) > 1) {
  combined_plot <- (
    (kinship_matrix / pop_meta / spp_meta / dendrogram) +
      plot_layout(guides = "collect", heights = c(37, 0.5, 0.5, 2)
      ) |
      (plot_spacer() / kinship_key / plot_spacer()) +
        plot_layout(heights = c(0.01, 36.6, 3.4)
        )
  ) + plot_layout(widths = c(39, 1))
}

combined_plot <- combined_plot &
  theme(
    legend.position = "left",
    legend.justification = "top",
    legend.key.size = unit(2, "mm"),
    legend.key.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(t = 2.1, b = -2.5, l = 0, r = -23, unit = "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold")
  )

ggsave(
  plot = combined_plot,
  filename = paste0(name, ".png"),
  dpi = 600,
  width  = 6.75 / 2,
  height = 6.75 / 2.1,
  bg = "white"
)
