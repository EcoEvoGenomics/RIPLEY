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

# Missingness is reported as a fraction, so it is rescaled for percentage axes
imiss$PCT_MISS <- 100 * imiss$F_MISS
lmiss$PCT_MISS <- 100 * lmiss$F_MISS
hwe$NEG_LOG10_P <- -log10(hwe$P_HWE)

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

# Ordered by filtering workflow, ending with F as an interpretive metric
panels <- list(
  list(data = idepth, column = "MEAN_DEPTH",
       title = "Sample Mean Sequencing Depth", xlab = ""),
  list(data = ldepth, column = "MEAN_DEPTH",
       title = "Site Mean Sequencing Depth (0th - 99th Percentile)", xlab = "",
       subset = ldepth$MEAN_DEPTH <= quantile(ldepth$MEAN_DEPTH, 0.99)),
  list(data = imiss, column = "PCT_MISS",
       title = "Sites Missing per Sample", xlab = "%"),
  list(data = lmiss, column = "PCT_MISS",
       title = "Samples Missing per Site", xlab = "%"),
  list(data = lqual, column = "QUAL",
       title = "Site Quality (0th - 99th Percentile)",
       xlab = expression(bolditalic("PHRED") ~ bold("Score")),
       subset = lqual$QUAL <= quantile(lqual$QUAL, 0.99)),
  list(data = frq, column = "MAF",
       title = "Site Minor Allele Frequency", xlab = ""),
  list(data = hwe, column = "NEG_LOG10_P",
       title = "Site Deviation from Hardy-Weinberg Equilibrium",
       xlab = expression(bolditalic(-log) * bold(""[10] ~ (P)))),
  list(data = het, column = "F",
       title = "Sample Inbreeding Coefficient", xlab = "F")
)

combined_plot <- panels |>
  lapply(\(panel) draw_histogram(
    panel$data,
    panel$column,
    panel$title,
    panel$xlab,
    subset = panel$subset
  )) |>
  wrap_plots(ncol = 2)

ggsave(
  plot = combined_plot,
  filename = paste0(name, ".png"),
  dpi = 600,
  width = 6.75,
  height = 6.75 / 2,
  bg = "white"
)
