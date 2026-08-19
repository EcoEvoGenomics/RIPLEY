library(tidyverse)
library(pheatmap)

args <- commandArgs(trailing = TRUE)
data <- read.table(args[1], header = TRUE)
name <- basename(args[1])

data <- data |>
  select(INDV1, INDV2, RELATEDNESS_PHI) |>
  pivot_wider(names_from = INDV2, values_from = RELATEDNESS_PHI) |>
  column_to_rownames("INDV1")

second_degree_relatedness <- 1 / (2^(7 / 2))
highlights <- as.matrix(data)
highlights <- formatC(highlights, format = "g", digits = 2)
highlights[is.na(highlights) | highlights <= second_degree_relatedness] <- ""

kinship_matrix <- pheatmap(
  mat = data,
  color = viridisLite::inferno(100),
  clustering_method = "ward.D2",
  border_color = NA,
  legend = FALSE,
  display_numbers = highlights,
  number_color = "red",
  fontsize_number = 5
)

n_inds <- nrow(data)
inches_per_ind <- 0.25

ggsave(
  plot = kinship_matrix,
  filename = paste(name, ".png", sep = ""),
  dpi = 600,
  height = n_inds * inches_per_ind,
  width = n_inds * inches_per_ind,
  bg = "white"
)
