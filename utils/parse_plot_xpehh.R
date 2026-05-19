# All hope abandon, ye who enter here ...

suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(forcats))
suppressMessages(library(gtools))

# -------- Parse command line arguments ----------------------------------------
args <- commandArgs(trailing = TRUE)
if (length(args) != 8) {
  stop("USAGE: Rscript plot_xpehh.R
  <path to list of xpehh.csv file paths>
  <path to list of xpehh.cand.csv file paths>
  <path to annotations .gff>
  <path to chromosome name conversion .tsv>
  <base p-value used to call candidate regions (e.g. 0.01)>
  <main plot width in mm>
  <candidate region plots width in mm>
  <plot height per scan in mm>"
  )
}
input_scans <- args[1]
input_cands <- args[2]
input_gff <- args[3]
chr_conversion_table <- args[4]
base_pval <- as.double(args[5])
width_mm <- as.integer(args[6])
cand_mm <- as.integer(args[7])
height_mm <- as.integer(args[8])

scan_files <- readLines(input_scans)
cand_files <- readLines(input_cands)

n_scans <- length(scan_files)
n_cands <- length(cand_files)
stopifnot(n_cands == n_scans)

# -------- Parse chromosomes rename tsv ----------------------------------------
chr_conversions <- read.table(chr_conversion_table)
renamed_chrs <- chr_conversions$V1
names(renamed_chrs) <- chr_conversions$V2

# -------- Parse input xp-EHH .CSV files ---------------------------------------
print("Parsing input files ...")

scan_list <- vector("list", length(scan_files))
cand_list <- vector("list", length(cand_files))

for (scan_index in seq_along(scan_files)) {

  file_path <- scan_files[scan_index]
  file_nrows <- as.integer(system(paste("wc -l <", file_path), intern = TRUE))
  scan_nsnps <- file_nrows - 1

  tmp <- paste(tempfile(), ".csv.tmp", sep = "")
  awk_fields <- "'{print $1,$2,$7,$8}'"
  paste("cat", file_path, "| awk -F,", awk_fields, ">", tmp) |> system()

  parsed_scanfile <- read.table(tmp, header = TRUE) |>
    as_tibble() |>
    rename(XPEHH = 3) |>
    drop_na(LOGPVALUE) |>
    mutate(SCAN = file_path) |>
    mutate(SCAN_BONFERRONI_THRESHOLD = -log10(base_pval / scan_nsnps))

  scan_list[[scan_index]] <- parsed_scanfile
  rm(parsed_scanfile)
  invisible(gc())

}
parsed_scans <- bind_rows(scan_list)

for (cand_index in seq_along(cand_files)) {

  file_path <- cand_files[cand_index]
  scan_index <- cand_index # Must be ordered identically in input
  matched_scan_path <- scan_files[scan_index]

  parsed_candfile <- read.csv(file_path, header = TRUE) |>
    as_tibble() |>
    mutate(SCAN = matched_scan_path) |>
    select(SCAN, CHR, START, END, MEAN_MRK) |>
    rename(MEAN_LOGPVALUE = MEAN_MRK) |>
    rowwise() |>
    mutate(
      MEAN_XPEHH = mean(
        parsed_scans$XPEHH[
          parsed_scans$SCAN == SCAN &
          parsed_scans$CHR == CHR &
          parsed_scans$POSITION >= START &
          parsed_scans$POSITION <= END
        ],
        na.rm = TRUE
      )
    ) |>
    ungroup()

  cand_list[[cand_index]] <- parsed_candfile
  rm(parsed_candfile)
  invisible(gc())

}
parsed_cands <- bind_rows(cand_list)

format_scans_cands <- function(scans_or_cands, renamed_chrs) {
  scans_or_cands <- scans_or_cands |>
    mutate(SCAN = factor(SCAN, levels = unique(SCAN))) |>
    mutate(CHR = factor(CHR, levels = gtools::mixedsort(unique(CHR)))) |>
    mutate(across(CHR, \(x) fct_recode(x, !!!renamed_chrs)))
}

parsed_scans <- format_scans_cands(parsed_scans, renamed_chrs)
parsed_cands <- format_scans_cands(parsed_cands, renamed_chrs)

write.csv(
  parsed_cands |>
    arrange(desc(abs(MEAN_XPEHH) * MEAN_LOGPVALUE)) |>
    mutate(SCAN = tools::file_path_sans_ext(basename(as.character(SCAN)))) |>
    mutate(SCAN = tools::file_path_sans_ext(SCAN)),
  row.names = FALSE,
  file = "ranked_cands.csv"
)
print("Input files parsed")

# -------- Parse annotations from GFF file -------------------------------------
print("Parsing annotations from GFF ...")

annots <- input_gff |>
  read_tsv(
    comment = "#",
    show_col_types = FALSE,
    col_names = c(
      "CHR", "SOURCE", "FEATURE",
      "START", "END", "SCORE", "STRAND", "FRAME",
      "ATTRIBUTE"
    )
  ) |>
  filter(FEATURE == "gene") |>
  mutate(CHR = factor(CHR)) |>
  mutate(across(CHR, \(x) fct_recode(x, !!!renamed_chrs)))

print("GFF parsed")

# -------- Plotting constants --------------------------------------------------
textsize <- 6
scan_names <- unique(parsed_scans$SCAN)
scan_labels <- tools::file_path_sans_ext(basename(scan_files))
scan_labels <- tools::file_path_sans_ext(scan_labels)
scan_labels_for_paths <- scan_labels
scan_labels <- str_replace_all(scan_labels, "_", " - ")
names(scan_labels_for_paths) <- scan_names
names(scan_labels) <- scan_names
chr_names <- levels(parsed_scans$CHR)
chr_labels_main <- chr_names
chr_labels_cand <- chr_names
for (chr_index in seq_along(chr_names)) {
  chr <- chr_labels_main[chr_index]
  chr_size <- max(parsed_scans$POSITION[parsed_scans$CHR == chr])
  if (nchar(chr) > (((chr_size / 1e9) * width_mm / 1.75))) {
    chr_labels_main[chr_index] <- ""
  }
  chr_labels_cand[chr_index] <- chr_conversions$V3[which(chr_conversions$V2 == as.character(chr))]
}
names(chr_labels_main) <- chr_names
names(chr_labels_cand) <- chr_names
facet_labels <- c(chr_labels_main, scan_labels)

compared_pops <- parsed_scans |>
  distinct(SCAN) |>
  mutate(SCAN_LABEL = scan_labels) |>
  group_by(SCAN) |>
  summarise(
    POP_A = stringr::str_split(as.character(SCAN_LABEL), " - ")[[1]][1],
    POP_B = stringr::str_split(as.character(SCAN_LABEL), " - ")[[1]][2],
    .groups = "drop"
  )

theme_common <- theme(
  axis.line = element_line(colour = "black", linewidth = 0.1),
  axis.text = element_text(colour = "black", size = textsize),
  axis.ticks.length = unit(0.5, "mm"),
  axis.ticks = element_line(colour = "black", linewidth = 0.1),
  axis.title.y = element_text(
    angle = 0,
    hjust = 0,
    vjust = 1,
    colour = "black",
    size = textsize,
    face = "bold"
  ),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.spacing.x = unit(0, "mm"),
  strip.background = element_blank(),
  strip.clip = "off",
  strip.text.x = element_text(
    colour = "black",
    size = textsize,
    margin = margin(0, 0, 0, 0, unit = "mm")
  ),
  strip.text.y = element_blank()
)

# -------- Plot main Manhattan plot --------------------------------------------
print("Plotting main plot ...")

manhattan_ymax <- ceiling(max(parsed_scans$LOGPVALUE))

chr_ranges <- parsed_scans |>
  group_by(CHR) |>
  summarise(
    xmin = min(POSITION),
    xmax = max(POSITION),
    is_last_chr = as.logical(unique(CHR == tail(levels(parsed_scans$CHR), n = 1))),
    .groups = "drop"
  )

manhattan <- ggplot(parsed_scans, aes(x = POSITION, y = LOGPVALUE)) +
  facet_grid(
    cols = vars(CHR),
    rows = vars(SCAN),
    labeller = as_labeller(facet_labels),
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  geom_rect(
    data = parsed_cands,
    aes(
      xmin = START,
      xmax = END,
      ymax = manhattan_ymax,
      ymin = 0
    ),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.1
  ) +
  geom_segment(
    data = chr_ranges,
    aes(
      x = xmax,
      xend = xmax,
      y = 0,
      yend = manhattan_ymax,
      alpha = is_last_chr
    ),
    colour = "black",
    linewidth = 0.1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = manhattan_ymax,
    colour = "black",
    linewidth = 0.1,
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.1,
  ) +
  geom_point(
    size = 0.25, stroke = 0, show.legend = FALSE
  ) +
  geom_point(
    aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
    colour = "red", size = 0.25, stroke = 0, show.legend = FALSE
  ) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(
    guide = guide_axis(cap = TRUE),
    breaks = seq(
      from = 0,
      to = manhattan_ymax,
      length.out = ceiling(height_mm / 7.5)
    ),
    labels = seq(
      from = 0,
      to = manhattan_ymax,
      length.out = ceiling(height_mm / 7.5)
    ) |> round(digits = 0),
    limits = c(0, manhattan_ymax)
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  ylab(expression(bold(bolditalic("-log")[10] ~ "P"))) +
  theme_common +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing.y = unit(height_mm / 7.5, "mm")
  )

cand_summary_ymax <- ceiling(max(parsed_cands$MEAN_XPEHH))
cand_summary_xmax <- ceiling(max(parsed_cands$MEAN_LOGPVALUE))
cand_summary_annotations <- tibble(
  SCAN = head(levels(parsed_cands$SCAN), n = 1),
  annot = expression(bold(mu * "(" * bolditalic("-log")[10] ~ "P)"))
)

cand_summary_plot <- parsed_cands |>
  ggplot(aes(y = MEAN_XPEHH, x = MEAN_LOGPVALUE)) +
  facet_grid(
    rows = vars(SCAN),
    labeller = as_labeller(facet_labels)
  ) +
  annotate(
    "segment",
    x = 0, xend = cand_summary_xmax,
    y = 0, yend = 0,
    colour = "black",
    linewidth = 0.1
  ) +
  geom_text(
    data = cand_summary_annotations,
    aes(label = annot),
    x = cand_summary_xmax,
    y = 0,
    vjust = 1.5,
    hjust = 1,
    size = 2
  ) +
  geom_text(
    data = compared_pops,
    aes(label = POP_A),
    x = cand_summary_xmax / 25,
    y = cand_summary_ymax,
    vjust = 1,
    hjust = 0,
    size = 2
  ) +
  geom_text(
    data = compared_pops,
    aes(label = POP_B),
    x = cand_summary_xmax / 25,
    y = -cand_summary_ymax,
    vjust = 0,
    hjust = 0,
    size = 2
  ) +
  geom_point(size = 0.5) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, cand_summary_xmax)
  ) +
  scale_y_continuous(
    guide = guide_axis(cap = TRUE),
    breaks = seq(
      from = -cand_summary_ymax,
      to = cand_summary_ymax,
      length.out = max(3, ceiling(height_mm / 7.5) + ((ceiling(height_mm / 7.5) %% 2) + 1))
    ),
    labels = seq(
      from = -cand_summary_ymax,
      to = cand_summary_ymax,
      length.out = max(3, ceiling(height_mm / 7.5) + ((ceiling(height_mm / 7.5) %% 2) + 1))
    ) |> round(digits = 0),
    limits = c(-cand_summary_ymax, cand_summary_ymax)
  ) +
  ylab(expression(bold(mu * "(" * bolditalic("xp") * "EHH)"))) +
  theme_common +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing.y = unit(height_mm / 7.5, "mm")
  )

outname <- paste(tools::file_path_sans_ext(input_scans), "_main.png", sep = "")
ggsave(
  (manhattan | cand_summary_plot) + plot_layout(widths = c(3, 1)),
  filename = outname,
  units = "mm",
  dpi = 1600,
  width = width_mm,
  height = height_mm * length(unique(parsed_scans$SCAN))
)

print(paste("Main plot saved to ", outname, sep = ""))

# -------- Plot candidate region Manhattan plots -------------------------------
print("Plotting candidate region plots and writing candidate genes to file ...")

candplot_padding <- 1e5 # How many basepairs to plot around candidate region
candplot_ymin <- -(manhattan_ymax / 6)
candplot_ygap <- candplot_ymin / 3

cand_outputdir <- tools::file_path_sans_ext(input_scans)
dir.create(cand_outputdir)

# First merge overlapping candidate intervals into combined plotting windows
merged_regions <- parsed_cands |>
  arrange(CHR, START, END) |>
  group_by(CHR) |>
  group_modify(\(cands, key) {

    merged <- list()
    current_start <- cands$START[1]
    current_end   <- cands$END[1]

    for (i in seq_len(nrow(cands))) {
      start <- cands$START[i]
      end <- cands$END[i]

      if (start <= current_end) {
        current_end <- max(current_end, end)
      } else {
        merged[[length(merged) + 1]] <- tibble(
          START = current_start,
          END = current_end
        )
        current_start <- start
        current_end <- end
      }
    }

    merged[[length(merged) + 1]] <- tibble(
      START = current_start,
      END = current_end
    )
    bind_rows(merged)

  }) |>
  ungroup() |>
  arrange(CHR, START, END)

for (region_index in seq_len(nrow(merged_regions))) {

  chr <- as.character(merged_regions$CHR[region_index])
  chr_min <- 0
  chr_max <- max(parsed_scans$POSITION[parsed_scans$CHR == chr])

  region_start <- merged_regions$START[region_index]
  region_end <- merged_regions$END[region_index]
  clamped_window_start <- max(chr_min, region_start - candplot_padding)
  clamped_window_end <- min(chr_max, region_end + candplot_padding)

  chr_outputdir <- paste(cand_outputdir, "/", chr, sep = "")
  dir.create(chr_outputdir)

  region_scans <- parsed_scans |>
    filter(CHR == chr) |>
    filter(POSITION >= clamped_window_start) |>
    filter(POSITION <= clamped_window_end) |>
    mutate(SCAN = factor(SCAN)) |>
    droplevels()

  region_cands <- parsed_cands |>
    filter(CHR == chr) |>
    filter(START >= clamped_window_start) |>
    mutate(END = ifelse(END > chr_max, chr_max, END)) |>
    filter(END <= clamped_window_end) |>
    mutate(SCAN = factor(SCAN, levels = levels(region_scans$SCAN))) |>
    droplevels()

  region_annots <- annots |>
    filter(CHR == chr) |>
    filter(START >= clamped_window_start) |>
    filter(END <= clamped_window_end)

  # Identify candidate genes and output them only if there are annotations
  if (nrow(region_annots) != 0) {

    region_annots <- region_annots |>
      crossing(SCAN = levels(region_scans$SCAN)) |>
      mutate(SCAN = factor(SCAN, levels = levels(region_scans$SCAN))) |>
      rowwise() |>
      mutate(
        IS_CANDGENE = {
          this_scan <- SCAN
          this_start <- START
          this_end <- END
          scan_cands <- region_cands |>
            filter(as.character(SCAN) == as.character(this_scan))
          no_scan_cands <- nrow(scan_cands) == 0
          if (no_scan_cands) {
            FALSE
          } else {
            any(
              this_start <= scan_cands$END &
                this_end >= scan_cands$START
            )
          }
        }
      ) |>
      ungroup()

    # Output gff file with candidate genes before modifying annotations for plot
    for (scan in levels(region_scans$SCAN)) {

      candgenes_in_scan <- region_annots |>
        filter(SCAN == scan, IS_CANDGENE) |>
        select(-SCAN, -IS_CANDGENE)

      no_candgenes <- nrow(candgenes_in_scan) == 0
      if (no_candgenes) {
        next
      }

      candgenes_outname <- paste(
        chr_outputdir, "/", chr,
        "_", format(region_start, scientific = FALSE),
        "_", format(region_end, scientific = FALSE),
        "_", scan_labels_for_paths[which(names(scan_labels_for_paths) == scan)],
        ".gff",
        sep = ""
      )

      candgenes_in_scan |> write_tsv(file = candgenes_outname)
    }

    # Reverse start and end for genes on minus strand for plotting arrows
    for (index in seq_len(nrow(region_annots))) {
      start <- region_annots$START[index]
      end <- region_annots$END[index]
      if (region_annots$STRAND[index] == "-") {
        region_annots$END[index] <- start
        region_annots$START[index] <- end
      }
    }

  }

  region_ymax <- ceiling(max(abs(region_scans$XPEHH)))
  region_ymax_padded <- region_ymax * 1.15
  region_ymin <- -region_ymax
  region_ymin_padded_axis <- -region_ymax_padded
  region_ymin_padded_limit <- region_ymin * 2
  region_yratio <- region_ymax / manhattan_ymax

  candhattan <- ggplot(region_scans, aes(x = POSITION, y = XPEHH)) +
    coord_cartesian(clip = "off") +
    facet_grid(
      rows = vars(SCAN),
      labeller = as_labeller(facet_labels)
    ) +
    annotate(
      "segment",
      x = clamped_window_start, xend = clamped_window_end,
      y = 0, yend = 0,
      colour = "black", linewidth = 0.1
    ) +
    geom_rect(
      data = region_cands,
      aes(
        xmin = START,
        xmax = END,
        ymax = region_ymax_padded,
        ymin = region_ymin_padded_axis
      ),
      inherit.aes = FALSE,
      fill = "red", colour = NA,
      alpha = 0.05
    ) +
    geom_text(
      data = compared_pops,
      aes(label = POP_A),
      x = clamped_window_start + (clamped_window_end - clamped_window_start) * 0.005,
      y = region_ymax_padded,
      vjust = 1,
      hjust = 0,
      size = 2
    ) +
    geom_text(
      data = compared_pops,
      aes(label = POP_B),
      x = clamped_window_start + (clamped_window_end - clamped_window_start) * 0.005,
      y = region_ymin_padded_axis,
      vjust = 0,
      hjust = 0,
      size = 2
    ) +
    geom_point(colour = "black", size = 0.75, stroke = 0, show.legend = FALSE) +
    geom_point(
      aes(
        alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD,
        size = LOGPVALUE
      ),
      colour = "red", stroke = 0,
      show.legend = FALSE
    ) +
    xlab(chr_labels_cand[which(names(chr_labels_cand) == chr)]) +
    scale_alpha_manual(values = c(0, 1)) +
    scale_size(range = c(0.5, 1.5)) +
    scale_x_continuous(
      breaks = seq(
        from = clamped_window_start,
        to = clamped_window_end,
        length.out = max(2, ceiling(width_mm / 40))
      ),
      labels = scales::label_number(
        accuracy = 0.1,
        scale    = 1 / 1e6,
        suffix   = " mbp"
      ),
      limits = c(
        clamped_window_start,
        clamped_window_end
      ),
      expand = expansion(mult = c(0, 0.05)),
      guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(
      guide = guide_axis(cap = TRUE),
      breaks = seq(
        from = region_ymin_padded_axis,
        to = region_ymax_padded,
        length.out = max(3, ceiling(height_mm / 7.5) + ((ceiling(height_mm / 7.5) %% 2) + 1))
      ),
      labels = seq(
        from = region_ymin_padded_axis,
        to = region_ymax_padded,
        length.out = max(3, ceiling(height_mm / 7.5) + ((ceiling(height_mm / 7.5) %% 2) + 1))
      ) |> round(digits = 0),
      limits = c(
        region_ymin_padded_limit + (candplot_ymin * region_yratio),
        region_ymax_padded
      )
    ) +
    ylab(expression(bold(bolditalic("xp") * "EHH"))) +
    theme_common +
    theme(
      axis.title.x = element_text(
        colour = "black",
        face = "bold",
        size = textsize,
        hjust = 0.48
      ),
      panel.spacing.y = unit(0, "mm"),
    )

  # Add annotation arrows if they exist in the window
  if (nrow(region_annots) != 0) {
    candhattan <- candhattan +
      geom_segment(
        data = region_annots,
        aes(
          x = START,
          xend = END,
          y = region_ymin_padded_limit + ((candplot_ymin - (candplot_ymin - candplot_ygap) / 2) * region_yratio),
          yend = region_ymin_padded_limit + ((candplot_ymin - (candplot_ymin - candplot_ygap) / 2) * region_yratio),
          colour = as.character(SCAN) == as.character(tail(levels(region_scans$SCAN), n = 1))
        ),
        arrow = arrow(
          angle = 30,
          length = unit(0.75, "mm"),
          type = "closed"
        ),
        linewidth = 0.25,
        lineend = "square",
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      scale_colour_manual(values = c(NA, "blue"))
  }

  candplot_outname <- paste(
    chr_outputdir, "/",
    chr, "_", format(region_start, scientific = FALSE),
    "_", format(region_end, scientific = FALSE), ".png",
    sep = ""
  )
  ggsave(
    candhattan,
    filename = candplot_outname,
    dpi = 1600,
    units = "mm",
    width = cand_mm,
    height = height_mm * length(unique(parsed_scans$SCAN))
  )
  print(
    paste(
      "Candidate region plot ",
      region_index, "/", nrow(merged_regions),
      " saved to ", candplot_outname,
      " and candidate genes extracted from your .gff",
      sep = ""
    )
  )

}

print("Done")
