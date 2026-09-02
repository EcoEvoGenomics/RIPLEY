library(tidyverse)
library(pheatmap)

args <- commandArgs(trailing = TRUE)
name <- basename(args[1])
data <- read.table(args[1], header = TRUE)
meta <- read.table(
  args[2],
  sep = ",",
  header = FALSE,
  col.names = c("ID", "Species", "Population", "Sex")
)

data <- data |>
  select(INDV1, INDV2, RELATEDNESS_PHI) |>
  pivot_wider(names_from = INDV2, values_from = RELATEDNESS_PHI) |>
  column_to_rownames("INDV1")

meta <- meta |>
  column_to_rownames("ID") |>
  select(Species, Population) |>
  select(where(\(x) n_distinct(x) > 1))

second_degree_relatedness <- 1 / (2^(7 / 2))
highlights <- as.matrix(data)
highlights <- formatC(highlights, format = "g", digits = 1)
highlights[is.na(highlights) | highlights <= second_degree_relatedness] <- ""

for (index in seq_len(nrow(highlights))) {
  highlights[index, index] <- ""
}

kinship_matrix <- pheatmap(
  mat = data,
  color = viridisLite::inferno(100),
  treeheight_col = 0,
  clustering_distance_rows = dist(1 - data),
  clustering_distance_cols = dist(1 - data),
  clustering_method = "ward.D2",
  annotation_row = meta,
  annotation_names_row = FALSE,
  display_numbers = highlights,
  number_color = "black",
  fontsize = 12,
  fontsize_number = 9,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend = FALSE
)

n_inds <- nrow(data)
inches_per_ind <- 0.25

ggsave(
  plot = kinship_matrix,
  filename = paste(name, ".png", sep = ""),
  dpi = 600,
  height = n_inds * 0.85 *  inches_per_ind,
  width = n_inds * inches_per_ind,
  bg = "white"
)
