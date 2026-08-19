library(ggplot2)

args <- commandArgs(trailing = TRUE)
data <- read.table(args[1], header = TRUE)
name <- basename(args[1])

ordered_chroms <- unique(data$CHR)
chrom_count <- length(ordered_chroms)
data$CHR <- factor(data$CHR, levels = ordered_chroms)

decay_plot <- ggplot(data, aes(x = DIST, y = AVG_R2)) +
  geom_point(alpha = 0.1, stroke = NA) +
  facet_wrap(vars(CHR)) +
  scale_y_continuous(
    name = bquote("Linkage diseqiulibrium, mean" ~ R^2),
    n.breaks = 3,
    limits = c(0, 1),
    expand = FALSE
  ) +
  scale_x_continuous(
    name = "Distance (bp)",
    n.breaks = 3
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(10, "mm")
  )

ggsave(
  plot = decay_plot,
  filename = paste(name, ".png", sep = ""),
  height = chrom_count / 2,
  width = chrom_count / 2,
  dpi = 600,
  bg = "white"
)
