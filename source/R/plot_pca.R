# Based on our PCA parsing scripts at https://github.com/EcoEvoGenomics/eegr

library(tidyverse)

args <- commandArgs(trailing = TRUE)

eigenvalues <- scan(args[1], quiet = TRUE)
variance_explained <- eigenvalues / sum(eigenvalues)
variance_explained <- signif(variance_explained, digits = 3)
variance_percent <- variance_explained * 100

eigenvectors <- read.table(args[2], header = FALSE)
eigenvectors_id_index  <- 1
eigenvectors_pid_index <- 2
eigenvectors <- eigenvectors[-eigenvectors_pid_index]
n_pc <- ncol(eigenvectors[-eigenvectors_id_index])
pc_indices <- seq(from = eigenvectors_id_index, to = n_pc) + 1
pc_names <- paste("PC", seq(n_pc), sep = "")
names(eigenvectors)[eigenvectors_id_index] <- "ID"
names(eigenvectors)[pc_indices] <- pc_names

meta <- read.table(
  args[3],
  sep = ",",
  header = FALSE,
  col.names = c("ID", "Species", "Population", "Sex")
)
n_species <- length(unique(meta$Species))
n_populations <- length(unique(meta$Population))

draw_pca <- function(eigenvectors, pcx_num, pcy_num, variance, meta, grouping) {

  pcx <- paste("PC", pcx_num, sep = "")
  pcy <- paste("PC", pcy_num, sep = "")
  pcx_variance <- variance[pcx_num]
  pcy_variance <- variance[pcy_num]

  expand <- 0.1
  xmax <- (round(max(eigenvectors[pcx_num + 1]) * 10) / 10) + expand
  xmin <- (round(min(eigenvectors[pcx_num + 1]) * 10) / 10) - expand
  ymax <- (round(max(eigenvectors[pcy_num + 1]) * 10) / 10) + expand
  ymin <- (round(min(eigenvectors[pcy_num + 1]) * 10) / 10) - expand

  text_size <- 2
  line_colour <- "black"
  point_size <- 2.5

  pca_plot <- eigenvectors |>
    left_join(meta, by = "ID") |>
    ggplot(
      aes(
        x = .data[[pcx]],
        y = .data[[pcy]],
        colour = .data[[grouping]],
        shape = Sex
      )
    ) +
    coord_equal(
      xlim = c(xmin, xmax), ylim = c(ymin, ymax),
      expand = FALSE, clip = "off"
    ) +
    geom_hline(yintercept = 0, colour = line_colour, linewidth = 0.1) +
    geom_vline(xintercept = 0, colour = line_colour, linewidth = 0.1) +
    annotate(
      "text", label = paste0(pcx, " (", pcx_variance, "%)"),
      x = xmax, y = (ymax - ymin) / 50,
      hjust = 1, vjust = 0,
      colour = line_colour, size = text_size
    ) +
    annotate(
      "text", label = xmax,
      x = xmax, y = (ymax - ymin) / -50,
      hjust = 1, vjust = 1,
      colour = line_colour, size = text_size
    ) +
    annotate(
      "text", label = xmin,
      x = xmin, y = (ymax - ymin) / -50,
      hjust = 0, vjust = 1,
      colour = line_colour, size = text_size
    ) +
    annotate(
      "text", label = paste0(pcy, " (", pcy_variance, "%)"),
      x = (xmax - xmin) / 50, y = ymax,
      hjust = 0, vjust = 1,
      colour = line_colour, size = text_size
    ) +
    annotate(
      "text", label = ymax,
      x = (xmax - xmin) / -50, y = ymax,
      hjust = 1, vjust = 1,
      colour = line_colour, size = text_size
    ) +
    annotate(
      "text", label = ymin,
      x = (xmax - xmin) / -50, y = ymin,
      hjust = 1, vjust = 0,
      colour = line_colour, size = text_size
    ) +
    geom_point(size = point_size) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.justification = "top",
      legend.key.spacing.y = unit(-2.5, "mm"),
      legend.margin = margin(t = 0, b = 2.5, unit = "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6, face = "bold"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank()
    )

  ggsave(
    plot = pca_plot,
    filename = paste(pcx, "_", pcy, "_", grouping, ".png", sep = ""),
    dpi = 600,
    width = 6.75,
    height = 6.75,
    bg = "white"
  )

}

plot_pcs <- seq(from = 1, to =  n_populations - 1)
for (pc in plot_pcs) {

  pcx <- pc
  pcy <- pcx + 1

  if (pcy > max(plot_pcs)) break
  draw_pca(eigenvectors, pcx, pcy, variance_percent, meta, "Population")

  if (n_species <= 1) next
  draw_pca(eigenvectors, pcx, pcy, variance_percent, meta, "Species")

}
