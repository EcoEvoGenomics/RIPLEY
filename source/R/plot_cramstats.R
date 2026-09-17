library(tidyverse)
library(patchwork)

args <- commandArgs(trailing = TRUE)
name <- tools::file_path_sans_ext(basename(args[1]))

stats <- read.table(
  args[1],
  header = FALSE,
  sep = "\t",
  quote = "",
  col.names = c("ID", "METRIC", "VALUE"),
  stringsAsFactors = FALSE
) |>
  filter(!is.na(suppressWarnings(as.numeric(VALUE)))) |>
  mutate(VALUE = as.numeric(VALUE)) |>
  pivot_wider(id_cols = ID, names_from = METRIC, values_from = VALUE)

# SN reports counts, not proportions, so rates are taken against the sample total
stats$PCT_MAPPED <- 100 * stats$`reads mapped` / stats$`sequences`
stats$PCT_PROPERLY_PAIRED <- 100 * stats$`reads properly paired` / stats$`sequences`
stats$PCT_MQ0 <- 100 * stats$`reads MQ0` / stats$`sequences`

draw_histogram <- function(data, x, title, xlab, bins = 30, subset = NULL) {

  if (!is.null(subset)) data <- data[subset, ]

  ggplot(data, aes(x = .data[[x]])) +
    geom_histogram(bins = bins, fill = "#40468A", colour = "black", linewidth = 0.15) +
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

# Panels fill the combined grid row-wise, two per row
panels <- list(
  list(column = "PCT_MAPPED",
       title = "Reads Mapped per Sample", xlab = "%"),
  list(column = "PCT_PROPERLY_PAIRED",
       title = "Reads Properly Paired per Sample", xlab = "%"),
  list(column = "PCT_MQ0",
       title = "Reads Mapped with Zero Quality per Sample", xlab = "%"),
  list(column = "error rate",
       title = "Mismatches per Mapped Base", xlab = ""),
  list(column = "average quality",
       title = "Sample Mean Base Quality",
       xlab = expression(bolditalic("PHRED") ~ bold("Score"))),
  list(column = "average length",
       title = "Sample Mean Read Length", xlab = "bp"),
  list(column = "insert size average",
       title = "Sample Mean Insert Size", xlab = "bp"),
  list(column = "insert size standard deviation",
       title = "Sample Insert Size Standard Deviation", xlab = "bp")
)

combined_plot <- panels |>
  lapply(\(panel) draw_histogram(stats, panel$column, panel$title, panel$xlab)) |>
  wrap_plots(ncol = 2)

ggsave(
  plot = combined_plot,
  filename = paste0(name, ".png"),
  dpi = 600,
  width = 6.75,
  height = 6.75 / 2,
  bg = "white"
)
