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

# Missingness is reported as a fraction, so it is rescaled for percentage axes
imiss$PCT_MISS <- 100 * imiss$F_MISS
lmiss$PCT_MISS <- 100 * lmiss$F_MISS

hwe$NEG_LOG10_P <- -log10(hwe$P_HWE)

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

# To-Do: Add .frq-plot later
panels <- list(
  list(slug = "het", data = het, column = "F",
       title = "Sample Inbreeding Coefficient", xlab = "F"),
  list(slug = "hwe", data = hwe, column = "NEG_LOG10_P",
       title = "Site Deviation from Hardy-Weinberg Equilibrium",
       xlab = expression(bolditalic(-log) * bold(""[10] ~ (P)))),
  list(slug = "idepth", data = idepth, column = "MEAN_DEPTH",
       title = "Sample Mean Sequencing Depth", xlab = ""),
  list(slug = "imiss", data = imiss, column = "PCT_MISS",
       title = "Sites Missing per Sample", xlab = "%"),
  list(slug = "ldepth", data = ldepth, column = "MEAN_DEPTH",
       title = "Site Mean Sequencing Depth (0th - 99th Percentile)", xlab = "",
       subset = ldepth$MEAN_DEPTH <= quantile(ldepth$MEAN_DEPTH, 0.99)),
  list(slug = "lmiss", data = lmiss, column = "PCT_MISS",
       title = "Samples Missing per Site", xlab = "%"),
  list(slug = "lqual", data = lqual, column = "QUAL",
       title = "Site Quality (0th - 99th Percentile)",
       xlab = expression(bolditalic("PHRED") ~ bold("Score")),
       subset = lqual$QUAL <= quantile(lqual$QUAL, 0.99))
)

file_key <- str_split(statfiles[1], "_")[[1]][1]

for (panel in panels) {

  ggsave(
    plot = draw_histogram(
      panel$data,
      panel$column,
      panel$title,
      panel$xlab,
      subset = panel$subset
    ),
    filename = paste0(file_key, "_popwise.", panel$slug, ".png"),
    dpi = 600,
    width = (6.75 / 2),
    height = (6.75 / 12) * length(unique(idepth$POP)),
    bg = "white"
  )

}
