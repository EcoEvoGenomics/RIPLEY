library(tidyverse)
library(patchwork)

args <- commandArgs(trailing = TRUE)
hihet <- read.table(args[1], header = TRUE)
k <- basename(args[1])
meta <- read.table(
  args[2],
  sep = ",",
  header = FALSE,
  col.names = c("ID", "Species", "Population", "Sex")
)

plot_data <- hihet |>
  left_join(meta, by = "ID") |>
  separate(COMPARISON, into = c("P1", "P2"), sep = "-")

triangle_data <- plot_data |>
  distinct(P1, P2) |>
  crossing(
    data.frame(
      x = c(0, 1, 0.5, 0),
      y = c(0, 0, 1, 0)
    )
  )

hihet_plot <- plot_data |>
  ggplot(
    aes(
      x = HI,
      y = HET,
      fill = Population
    )
  ) +
  xlab("Hybrid Index") +
  ylab("Heterozygosity") +
  coord_equal() +
  facet_grid(cols = vars(P1), rows = vars(P2), switch = "x") +
  geom_polygon(
    data = triangle_data,
    aes(x = x, y = y),
    colour = "black",
    fill = NA,
    linewidth = 0.1,
    inherit.aes = FALSE
  ) +
  stat_function(
    fun = \(x) x * 2 * (1 - x),
    xlim = c(0, 1),
    linetype = 2,
    linewidth = 0.15,
    colour = "black"
  ) +
  geom_text(
    data = plot_data |> distinct(P1, P2) |> mutate(LABEL = P1),
    aes(
      label = LABEL,
      colour = factor(P1, levels = paste("p", seq(9), sep = "")) # Cols in palette
    ),
    x = 0, y = -0.1,
    hjust = 0.5,
    vjust = 0.5,
    size = 6 / .pt,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = plot_data |> distinct(P1, P2) |> mutate(LABEL = P2),
    aes(
      label = LABEL,
      colour = factor(P2, levels = paste("p", seq(9), sep = "")) # Cols in palette
    ),
    x = 1, y = -0.1,
    hjust = 0.5,
    vjust = 0.5,
    size = 6 / .pt,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(
    size = min(2.5, 2.5 / (0.75 * length(unique(plot_data$P1)))),
    pch = 21,
    stroke = 0.15
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 2.5)
    )
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(
    expand = expansion(add = 0.2),
    guide = guide_axis(cap = "both"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  scale_y_continuous(
    expand = expansion(add = 0.2),
    guide = guide_axis(cap = "both"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.15),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(size = 6, face = "bold"),
    axis.ticks = element_line(linewidth = 0.15),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification.inside = c(1, 1),
    legend.key.height = unit(2, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.margin = margin(t = 0, b = 0, unit = "mm"),
    legend.text = element_text(size = 6, hjust = 0),
    legend.title = element_text(size = 6, face = "bold", hjust = 0),
    strip.placement = "inside",
    strip.background.x = element_rect(fill = NA, colour = NA),
    strip.text.x = element_text(colour = NA, size = 0, margin = margin(t = 2, unit = "mm")),
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines")
  )

ggsave(
  plot = hihet_plot,
  filename = paste("k", k, "_hihet.png", sep = ""),
  dpi = 600,
  width = 6.75 / 2,
  height = 6.75 / 2
)
