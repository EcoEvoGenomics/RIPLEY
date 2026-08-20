library(tidyverse)

args <- commandArgs(trailing = TRUE)
chroms <- strsplit(args[2], ",")[[1]]
data <- read.table(args[1], header = TRUE) |> filter(CHROM %in% chroms)
chroms <- chroms[chroms %in% data$CHROM] # If expected chroms are not in data
name <- basename(args[1])
chrom_conversions <- read.table(args[3]) |> filter(V1 %in% chroms)

renamed_chroms <- chrom_conversions$V1
names(renamed_chroms) <- chrom_conversions$V2

n_chroms <- length(chroms)
bin_size <- data$BIN_START[2] - data$BIN_START[1]
xmax <- max(data$BIN_START) + bin_size

density_plot <- data |>
  ggplot(
    aes(
      x = BIN_START,
      y = factor(CHROM, levels = rev(chroms)),
      fill = log10(SNP_COUNT)
    )
  ) +
  geom_tile(height = 0.75) +
  scale_fill_viridis_c(
    name = expression(log[10] ~ (SNPs)),
    option = "magma"
  ) +
  scale_x_continuous(
    guide = guide_axis(cap = TRUE),
    expand = expansion(add = 0),
    limits = c((xmax * -0.015), (xmax * 1.03)),
    breaks = seq(from = 0, to = xmax, length.out = 3),
    position = "top",
    labels = scales::label_number(
      accuracy = 1,
      scale  = 1 / 1e6,
      suffix   = " mbp"
    )
  ) +
  scale_y_discrete(
    expand = 0.4 / n_chroms,
    labels = rev(names(renamed_chroms))
  ) +
  theme_void() +
  theme(
    axis.line.x = element_line(colour = "black", linewidth = 0.1),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.ticks.length = unit(0.5, "mm"),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
    legend.position = "bottom",
    legend.justification.bottom = c("right", "top"),
    legend.margin = margin(r = 0.5, unit = "cm"),
    legend.text = element_text(size = 6),
    legend.ticks = element_blank(),
    legend.title = element_text(size = 6, face = "bold", vjust = 0.85),
    legend.key.width  = unit(0.75, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

inches_per_chrom <- 0.125

ggsave(
  plot = density_plot,
  filename = paste(name, ".png", sep = ""),
  dpi = 600,
  height = (n_chroms * inches_per_chrom) + 1,
  width = 6.75,
  bg = "white"
)
