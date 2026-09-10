library(tidyverse)
library(patchwork)

args <- commandArgs(trailing = TRUE)
chroms <- strsplit(args[2], ",")[[1]]
data <- read.table(args[1], header = TRUE) |> filter(CHROM %in% chroms)
chroms <- chroms[chroms %in% data$CHROM] # If expected chroms are not in data
name <- basename(args[1])
chrom_labels <- read.table(args[3], sep = ",") |> filter(V1 %in% chroms)

renamed_chroms <- chrom_labels$V1
names(renamed_chroms) <- chrom_labels$V2

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
  ggtitle("Chromosome") +
  geom_tile(height = 0.75, colour = "black", fill = "black") +
  geom_tile(height = 0.75) +
  scale_fill_viridis_c(
    limits = c(0, max(log10(data$SNP_COUNT))),
    name = expression(log[10] ~ (SNPs)),
    option = "magma"
  ) +
  scale_x_continuous(
    guide = guide_axis(cap = TRUE),
    expand = expansion(add = 0),
    limits = c((xmax * -0.015), (xmax * 1.03)),
    breaks = seq(from = 0, to = xmax, length.out = 3),
    # position = "top",
    labels = scales::label_number(
      accuracy = 1,
      scale  = 1 / 1e6,
      suffix   = " mbp"
    )
  ) +
  scale_y_discrete(
    expand = expansion(add = 1),
    labels = rev(names(renamed_chroms))
  ) +
  theme_void() +
  theme(
    axis.line.x = element_line(colour = "black", linewidth = 0.1),
    axis.text.x = element_text(size = 5, margin = margin(t = 1, unit = "mm")),
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.ticks.length = unit(0.5, "mm"),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
    plot.title = element_text(size = 6, face = "bold", hjust = 0),
    legend.position = "none"
  )

density_key <- data |>
  ggplot(
    aes(
      y = 0,
      x = seq(
        0,
        max(log10(data$SNP_COUNT), na.rm = TRUE),
        length.out = length(data$SNP_COUNT)
      ),
      fill = seq(
        0,
        max(log10(data$SNP_COUNT), na.rm = TRUE),
        length.out = length(data$SNP_COUNT)
      )
    )
  ) +
  ggtitle(expression(bold("Number of SNPs (" ~ italic("log"[10]) ~ ")"))) +
  coord_cartesian(expand = FALSE) +
  geom_tile(show.legend = FALSE, height = 1.075, colour = "black") +
  geom_raster(show.legend = FALSE) +
  scale_x_continuous(
    expand = expansion(add = 0),
    limit = c(
      max(log10(data$SNP_COUNT), na.rm = TRUE) * - 0.015,
      max(log10(data$SNP_COUNT), na.rm = TRUE) * 1.03
    ),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5, margin = margin(t = 1, unit = "mm")),
    axis.ticks.length = unit(0.25, units = "mm"),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.15),
    plot.title = element_text(size = 6, face = "bold", hjust = 0, vjust = 4)
  )

inches_per_chrom <- 0.125

ggsave(
  plot = (density_key / plot_spacer() / density_plot) + plot_layout(heights = c(0.75, 0.75, n_chroms)),
  filename = paste(name, ".png", sep = ""),
  dpi = 600,
  height = (n_chroms * inches_per_chrom) + 1,
  width = 6.75,
  bg = "white"
)
