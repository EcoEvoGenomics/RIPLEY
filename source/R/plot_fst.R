library(tidyverse)

args <- commandArgs(trailing = TRUE)
fst <- read.table(
  args[1], header = FALSE,
  col.names = c("POP_A", "POP_B", "FST")
)

fst_symmetric <- fst |>
  bind_rows(rename(fst, POP_A = POP_B, POP_B = POP_A)) |>
  complete(POP_A, POP_B, fill = list(FST = 0))

distance <- fst_symmetric |>
  pivot_wider(names_from = POP_B, values_from = FST) |>
  column_to_rownames("POP_A") |>
  as.matrix()

distance <- distance[rownames(distance), rownames(distance)]
distance <- pmax(distance, 0) # Redundant in RUN_PAIRWISE_FST.nf
diag(distance) <- 0

clust <- hclust(as.dist(distance), method = "ward.D2")
pop_order <- clust$labels[clust$order]

fst <- fst_symmetric |>
  ggplot(
    aes(
      x = factor(POP_A, levels = pop_order),
      y = factor(POP_B, levels = pop_order),
      fill = FST,
      label = signif(FST, digits = 2)
    )
  ) +
  ggtitle(expression(bold("Weighted Mean F"[ST]))) +
  coord_equal(expand = FALSE) +
  geom_tile(show.legend = FALSE, colour = "black") +
  geom_text(size = 2) +
  scale_fill_distiller(palette = "Reds", direction = 2) +
  theme_void() +
  theme(
    axis.text.y = element_text(
      margin = margin(2.5, 2.5, 2.5, 0, unit = "mm"),
      size = 9,
      hjust = 1
    ),
    axis.text.x = element_text(
      margin = margin(2.5, 2.5, 0, 2.5, unit = "mm"),
      size = 9,
      angle = 90,
      hjust = 1
    ),
    panel.grid.major.x = element_line(colour = "grey75", linetype = 2),
    plot.title = element_text(
      margin = margin(0, 0, 2.5, 0, unit = "mm")
    )
  )

fullwidth_inches <- 6.75
plotwidth_inches <- ifelse(
  nrow(distance) > 6,
  fullwidth_inches,
  fullwidth_inches / 2
)

ggsave(
  plot = fst,
  filename = "weighted_mean_fst.png",
  bg = "white",
  dpi = 600,
  width = plotwidth_inches,
  height = plotwidth_inches,
  units = "in",
)
