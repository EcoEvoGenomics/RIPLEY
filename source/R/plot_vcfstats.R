library(tidyverse)
library(patchwork)

args <- commandArgs(trailing = TRUE)
name <- tools::file_path_sans_ext(basename(args[1]))
frq_path <- args[1]
idepth <- read.table(args[2], header = TRUE)
imiss <- read.table(args[3], header = TRUE)
ldepth <- read.table(args[4], header = TRUE)
lqual <- read.table(args[5], header = TRUE)
lmiss <- read.table(args[6], header = TRUE)
het <- read.table(args[7], header = TRUE)
hwe <- read.table(args[8], header = TRUE)

# Count greatest allele number to format frq table for MAF
frq_header_index <- 1
frq_field_counts <- count.fields(frq_path)[-frq_header_index]
frq_n_metadata_columns <- 4
frq_max_alleles <- max(frq_field_counts, na.rm = TRUE) - frq_n_metadata_columns
frq_col_names <- c(
  "CHROM",
  "POS",
  "N_ALLELES",
  "N_CHROMS",
  paste0("A", seq_len(frq_max_alleles))
)

frq <- read.table(
  frq_path,
  skip = frq_header_index,
  header = FALSE,
  fill = TRUE,
  col.names = frq_col_names,
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

frq$MAF <- frq[grep("^A", names(frq), value = TRUE)] |> apply(1, \(x) min(x))

draw_histogram <- function(data, x, title, xlab, bins = 30, subset = NULL) {

  if (!is.null(subset)) data <- data[subset, ]

  ggplot(data, aes(x = .data[[x]])) +
    geom_histogram(bins = bins, fill = "#721494", colour = "black", linewidth = 0.15) +
    scale_x_continuous(
      breaks = \(lim) pretty(lim, n = bins / 3 |> round()),
      expand = expansion(mult = c(0.025, 0)),
      guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(
      guide = guide_axis(cap = "both"),
      breaks = \(lim) pretty(lim, n = 2) |> round() |> unique(),
      labels = scales::label_number(accuracy = 1)
    ) +
    ggtitle(title) +
    xlab(xlab) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 6,
        face = "bold",
        hjust = 0,
        margin = margin(t = -1.25, b = 1, unit = "mm")
      ),
      axis.text = element_text(size = 5, colour = "black"),
      axis.title.x = element_text(size = 5, colour = "black", face = "bold", margin = margin(t = 1, b = 3, unit = "mm")),
      axis.title.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(colour = "black", linewidth = 0.15),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 0.15),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )

}

p1 <- draw_histogram(
  het, "F",
  "Sample Inbreeding Coefficient",
  "F"
)

p2 <- draw_histogram(
  idepth, "MEAN_DEPTH",
  "Sample Mean Sequencing Depth",
  ""
)

p3 <- draw_histogram(
  imiss,
  "F_MISS",
  "Fraction of Sites Missing per Sample",
  ""
)

p4 <- draw_histogram(
  ldepth, "MEAN_DEPTH",
  "Site Mean Sequencing Depth (0th - 99th Percentile)",
  "",
  subset = ldepth$MEAN_DEPTH <= quantile(ldepth$MEAN_DEPTH, 0.99)
)

p5 <- draw_histogram(
  lqual, "QUAL",
  "Site Quality (0th - 99th Percentile)",
  expression(bolditalic("PHRED") ~ bold("Score")),
  subset = lqual$QUAL <= quantile(lqual$QUAL, 0.99)
)

p6 <- draw_histogram(
  lmiss, "F_MISS",
  "Fraction of Samples Missing per Site",
  ""
)

p7 <- draw_histogram(
  frq, "MAF",
  "Site Minor Allele Frequency",
  ""
)

hwe$NEG_LOG10_P <- -log10(hwe$P_HWE)
p8 <- draw_histogram(
  hwe, "NEG_LOG10_P",
  "Site Deviation from Hardy-Weinberg Equilibrium",
  expression(bolditalic(-log) * bold(""[10] ~ (P)))
)

combined_plot <- (p2 | p4) / (p3 | p6) / (p1 | p5) / (p7 | p8)

ggsave(
  plot = combined_plot,
  filename = paste0(name, ".png"),
  dpi = 600,
  width = 6.75,
  height = 6.75 / 2,
  bg = "white"
)
