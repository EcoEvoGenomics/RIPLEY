library(tidyverse)
library(patchwork)
library(tools)

args <- commandArgs(trailing = TRUE)
statdir <- args[1]
statfiles <- dir(statdir)
het <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "het")]),
  header = TRUE
)
hwe <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "hwe")]),
  header = TRUE
)
idepth <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "idepth")]),
  header = TRUE
)
imiss <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "imiss")]),
  header = TRUE
)
ldepth <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "ldepth")]),
  header = TRUE
)
lmiss <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "lmiss")]),
  header = TRUE
)
lqual <- read.table(
  paste0(statdir, "/", statfiles[which(file_ext(statfiles) == "lqual")]),
  header = TRUE
)

draw_histogram <- function(data, x, title, xlab, bins = 30, subset = NULL) {

  library(ggplot2)

  if (!is.null(subset)) data <- data[subset, ]

  ggplot(data, aes(x = .data[[x]])) +
    facet_grid(
      rows = vars(.data[["POP"]])
    ) +
    geom_histogram(
      bins = bins,
      fill = "#721494",
      colour = "black",
      linewidth = 0.15
    ) +
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
      axis.title.x = element_text(
        size = 5,
        colour = "black",
        face = "bold",
        margin = margin(t = 1, b = 3, unit = "mm")
      ),
      axis.title.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(colour = "black", linewidth = 0.15),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 0.15),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      strip.clip = "off",
      strip.background = element_blank(),
      strip.text.y.right = element_text(
        size = 5,
        colour = "black",
        hjust = 0,
        vjust = -0.1,
        angle = 0
      )
    )

}

# p1 would be for .frq

p2 <- draw_histogram(
  het, "F",
  "Sample Inbreeding Coefficient",
  "F"
)

hwe$NEG_LOG10_P <- -log10(hwe$P_HWE)
p3 <- draw_histogram(
  hwe, "NEG_LOG10_P",
  "Site Deviation from Hardy-Weinberg Equilibrium",
  expression(bolditalic(-log) * bold(""[10] ~ (P)))
)

p4 <- draw_histogram(
  idepth, "MEAN_DEPTH",
  "Sample Mean Sequencing Depth",
  ""
)

p5 <- draw_histogram(
  imiss,
  "F_MISS",
  "Fraction of Sites Missing per Sample",
  ""
)

p6 <- draw_histogram(
  ldepth, "MEAN_DEPTH",
  "Site Mean Sequencing Depth (0th - 99th Percentile)",
  "",
  subset = ldepth$MEAN_DEPTH <= quantile(ldepth$MEAN_DEPTH, 0.99)
)

p7 <- draw_histogram(
  lmiss, "F_MISS",
  "Fraction of Samples Missing per Site",
  ""
)

p8 <- draw_histogram(
  lqual, "QUAL",
  "Site Quality (0th - 99th Percentile)",
  expression(bolditalic("PHRED") ~ bold("Score")),
  subset = lqual$QUAL <= quantile(lqual$QUAL, 0.99)
)

for (idx in seq_along(statfiles)) {

  file_key <- str_split(statfiles[idx], "_")[[1]][1]
  file_ext <- file_ext(statfiles[idx])

  if (file_ext == "frq") next # To-Do: Add .frq-plot later

  ggsave(
    plot = get(paste0("p", idx)),
    filename = paste0(file_key, "_popwise.", file_ext, ".png"),
    dpi = 600,
    width = (6.75 / 2),
    height = (6.75 / 12) * length(unique(idepth$POP)),
    bg = "white"
  )
}
