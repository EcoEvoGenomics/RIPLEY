library(tidyverse)
library(patchwork)

args <- commandArgs(trailing = TRUE)
dcms_dir <- args[1L]

# Only the first two columns are named, so further columns may be added later
read_group_metadata <- function(path) {
    group_metadata <- read.table(
        path,
        sep = ",",
        header = FALSE,
        comment.char = "" # Avoids hex codes reading in as comments
    )
    names(group_metadata)[1L:2L] <- c("GROUP", "COLOUR")
    return(group_metadata)
}

population_metadata <- read_group_metadata(args[2L])
chrom_labels <- read.table(args[3L], sep = ",", header = FALSE)
names(chrom_labels)[1L:2L] <- c("CHROM", "LABEL")

# DCMS is published one file per population pair, each spanning the whole genome
dcms_paths <- list.files(dcms_dir, pattern = "\\.dcms\\.tsv$", full.names = TRUE)
dcms <- map_dfr(
    dcms_paths,
    function(path) {
        populations <- str_remove(basename(path), "\\.dcms\\.tsv$") |> str_split_1("_")
        table <- read.delim(path, na.strings = c("NA", "nan", "")) |>
            mutate(POP_A = populations[1L], POP_B = populations[2L])
        return(table)
    }
)

dcms <- dcms |>
    mutate(
        # Starts are 1-based, so the half-open centre keeps the raster tiles inside the axis
        MIDPOINT = (WINDOW_START - 1L + WINDOW_END) / 2L,
        CHR = factor(CHR, levels = chrom_labels$CHROM),
        POP_A = factor(POP_A, levels = population_metadata$GROUP),
        POP_B = factor(POP_B, levels = population_metadata$GROUP),
        COMPARISON = fct_inorder(paste(POP_A, POP_B, sep = " - "))
    )

# Contigs are drawn end to end on one continuous axis, so each needs an offset
genome <- dcms |>
    summarise(LENGTH = max(WINDOW_END), .by = CHR) |>
    arrange(CHR) |>
    mutate(
        OFFSET = lag(cumsum(LENGTH), default = 0L),
        CENTRE = OFFSET + LENGTH / 2L,
        LABEL = chrom_labels$LABEL[match(as.character(CHR), chrom_labels$CHROM)],
        SHADE = row_number() %% 2L == 0L # Alternating tone stands in for contig boundaries
    )

dcms <- dcms |>
    left_join(select(genome, CHR, OFFSET, SHADE), by = "CHR") |>
    mutate(GENOME_POSITION = MIDPOINT + OFFSET)

shared_keys <- c(
    "CHR",
    "WINDOW_START",
    "WINDOW_END",
    "POP_A",
    "POP_B",
    "MIDPOINT",
    "COMPARISON",
    "OFFSET",
    "SHADE",
    "GENOME_POSITION"
)

# Available statistics vary with which scans ran
stats <- setdiff(names(dcms), c(shared_keys, "DCMS"))
stats <- stats[!str_ends(stats, "_P")]

# Presentation is hardcoded so a newly added statistic is labelled deliberately;
# anything unlisted falls back to its raw column name.
stat_catalogue <- list(
    pi = list(title = "Nucleotide diversity", symbol = quote(pi), short = quote(pi)),
    watterson_theta = list(
        title = "Watterson's estimator",
        symbol = quote(theta[W]),
        short = quote(theta[W])
    ),
    tajima_d = list(title = "Tajima's", symbol = quote(D), short = quote(D)),
    ihs = list(
        title = "Integrated haplotype score",
        symbol = quote(max ~ group("|", iHS, "|")),
        short = quote(iHS)
    ),
    dxy = list(title = "Absolute divergence", symbol = quote(d[XY]), short = quote(d[XY])),
    fst = list(title = "Fixation index", symbol = quote(F[ST]), short = quote(F[ST])),
    xpehh = list(
        title = "Extended haplotype homozygosity",
        symbol = quote(max ~ group("|", xpEHH, "|")),
        short = quote(xpEHH)
    )
)

track_title <- function(text, symbol) {
    if (is.null(symbol)) {
        return(bquote(bold(.(text))))
    }
    return(bquote(bold(.(text)) ~ italic(.(symbol))))
}

# calculate_dcms.R suffixes a statistic per population of the pair; the rest are pairwise
stat_role <- function(column) {
    role <- rep("PAIR", length(column))
    role[str_ends(column, "_A")] <- "A"
    role[str_ends(column, "_B")] <- "B"
    return(role)
}

stat_base <- function(column) {
    return(str_remove(column, "_[AB]$"))
}

# Rows are too short for a full title, so the conventional symbol stands alone
short_label <- function(column) {
    entry <- stat_catalogue[[stat_base(column)]]
    if (is.null(entry)) {
        return(bquote(.(column)))
    }
    return(bquote(italic(.(entry$short))))
}

# ---- Panel builders --------------------------------------------------------

# Row count depends on which scans ran, so the heatmap saturates instead of growing without bound.
# The constant is the fixed overhead of titles, the axis and the legends, so the two
# together cap the figure at 2.75 in however many statistics are supplied.
figure_height <- function(n_rows) {
    return(1.2 + 1.55 * (1.0 - exp(-0.18 * n_rows)))
}

# Loose breaks bracket the data, so pinning the limits to them caps every axis identically.
# They are resolved here rather than by the scale, which would re-derive them from the
# limits it was given and could widen the set beyond those limits.
axis_breaks <- function(values, n = 4L) {
    return(scales::extended_breaks(n = n, only.loose = TRUE)(range(values, na.rm = TRUE)))
}

# The lighter tone is the track colour blended halfway to white
shade_pair <- function(colour) {
    return(c(`FALSE` = colour, `TRUE` = colorRampPalette(c(colour, "white"))(4L)[2L]))
}

genome_scale <- function(show_axis) {
    return(scale_x_continuous(
        expand = expansion(add = 0L),
        limits = c(0L, max(genome$OFFSET + genome$LENGTH)),
        breaks = genome$CENTRE,
        labels = if (show_axis) genome$LABEL else NULL
    ))
}

base_theme <- function() {
    return(theme_classic(base_size = 6L) +
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_line(colour = "black", linewidth = 0.15),
            axis.ticks = element_line(colour = "black", linewidth = 0.15),
            axis.ticks.length = unit(0.5, "mm"),
            axis.text = element_text(size = 5L),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 6L),
            plot.title = element_text(size = 6L, hjust = 0L, margin = margin(b = 4L))
        ))
}

track <- function(data, column, colour, title) {
    breaks <- axis_breaks(data[[column]])
    plot <- ggplot(data, aes(x = GENOME_POSITION, y = .data[[column]])) +
        geom_line(
            aes(colour = SHADE, group = CHR), # Grouping also stops the line bridging contigs
            linewidth = 0.2
        ) +
        scale_colour_manual(values = shade_pair(colour), guide = "none") +
        ggtitle(title) +
        genome_scale(show_axis = TRUE) +
        scale_y_continuous(
            guide = guide_axis(cap = "both"),
            expand = expansion(add = 0L),
            limits = range(breaks),
            breaks = breaks
        ) +
        labs(x = "Chromosome", y = NULL) +
        base_theme() +
        theme(
            axis.title.x = element_text(margin = margin(t = 2L)),
            plot.margin = margin(b = 10L, r = 3L)
        )
    return(plot)
}

# One strip per statistic, sharing the axis drawn under the track above
heatmap_panel <- function(data, order) {
    plot <- ggplot(data, aes(x = GENOME_POSITION, y = STAT, fill = SIGNIFICANCE)) +
        geom_raster() +
        scale_fill_viridis_c(
            option = "inferno",
            direction = 1L,
            guide = guide_colourbar(
                title = expression(-log["10"] ~ italic(p)),
                barwidth = unit(16L, "mm"),
                barheight = unit(1.5, "mm"),
                ticks.linewidth = 0.3,
                frame.colour = "black",
                frame.linewidth = 0.15
            )
        ) +
        ggtitle(bquote(bold("Fractional rank") ~ bolditalic(p) ~ bold("values"))) +
        genome_scale(show_axis = FALSE) +
        scale_y_discrete(limits = order, expand = expansion(add = 0L)) +
        labs(x = NULL, y = NULL) +
        base_theme() +
        theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text = element_blank(),
            legend.title = element_text(size = 5L),
            legend.text = element_text(size = 5L),
            legend.position = "bottom",
            legend.justification = "left",
            legend.margin = margin(t = 2L),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.15),
            plot.margin = margin(r = 3L)
        )
    return(plot)
}

# Clustering interleaves the roles, so population is carried per row rather than by a spanning bar.
# The strip also carries the row labels, placing them between itself and the heatmap.
population_strip <- function(order, colours) {
    strip <- data.frame(STAT = order, ROLE = stat_role(order))
    row_labels <- map(order, short_label)
    names(row_labels) <- order
    plot <- ggplot(strip, aes(x = 0L, y = STAT, fill = ROLE)) +
        geom_tile() +
        scale_fill_manual(values = colours, guide = "none") +
        scale_y_discrete(
            limits = order,
            labels = function(breaks) {
                return(row_labels[breaks])
            },
            expand = expansion(add = 0L)
        ) +
        scale_x_continuous(expand = expansion(add = 0L)) +
        theme_void() +
        theme(
            axis.text.y = element_text(size = 5L, hjust = 1L, margin = margin(r = 3L)),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.15),
            plot.margin = margin(l = 3L)
        )
    return(plot)
}

# Axes are transposed so the tree opens leftward against the heatmap rows.
# Being the rightmost panel, it also carries the population key so that key ends flush with the
# figure edge; the tiles that declare it are given no area and so never draw.
dendrogram_panel <- function(clustering, colours, labels) {
    segments <- ggdendro::dendro_data(as.dendrogram(clustering))
    leaves <- mutate(segments$labels, ROLE = stat_role(as.character(label)))
    plot <- ggplot(ggdendro::segment(segments)) +
        geom_segment(aes(x = y, y = x, xend = yend, yend = xend), linewidth = 0.15) +
        geom_tile(
            data = leaves,
            mapping = aes(x = 0L, y = x, fill = ROLE),
            width = 0L,
            height = 0L
        ) +
        scale_fill_manual(
            values = colours,
            labels = labels,
            breaks = intersect(c("A", "B", "PAIR"), leaves$ROLE),
            guide = guide_legend(
                title = NULL,
                keywidth = unit(2L, "mm"),
                keyheight = unit(2L, "mm"),
                override.aes = list(colour = "black", linewidth = 0.15)
            )
        ) +
        scale_x_continuous(expand = expansion(add = 0L)) +
        scale_y_continuous(
            limits = c(0.5, nrow(segments$labels) + 0.5),
            expand = expansion(add = 0L)
        ) +
        theme_void() +
        theme(
            legend.text = element_text(size = 5L),
            legend.position = "bottom",
            legend.justification = "right",
            legend.margin = margin(t = 2L),
            plot.margin = margin()
        )
    return(plot)
}

# ---- Figure builder --------------------------------------------------------

draw_figure <- function(comparison) {
    data <- filter(dcms, COMPARISON == comparison)
    populations <- c(
        A = as.character(data$POP_A[1L]),
        B = as.character(data$POP_B[1L])
    )
    colours <- c(
        population_metadata$COLOUR[match(populations, population_metadata$GROUP)],
        PAIR = "grey30"
    )
    names(colours)[1L:2L] <- c("A", "B")
    labels <- c(
        A = populations[["A"]],
        B = populations[["B"]],
        PAIR = paste(populations[["A"]], "\u00d7", populations[["B"]])
    )

    # calculate_dcms.R publishes the rank-based p-value behind each statistic alongside it
    significance <- data[paste0(stats, "_P")] |>
        set_names(stats) |>
        mutate(across(everything(), \(p) -log10(p)))

    # Spearman mirrors the correlation DCMS itself weights by
    correlation <- cor(significance, method = "spearman", use = "pairwise.complete.obs")
    correlation[is.na(correlation)] <- 0L
    clustering <- hclust(as.dist(1L - correlation), method = "ward.D2")
    order <- clustering$labels[clustering$order]

    tiles <- significance |>
        mutate(GENOME_POSITION = data$GENOME_POSITION) |>
        pivot_longer(all_of(stats), names_to = "STAT", values_to = "SIGNIFICANCE")

    figure <- wrap_plots(
        track(
            data = data,
            column = "DCMS",
            colour = "black",
            title = track_title("Decorrelated composite signal", quote(DCMS))
        ),
        heatmap_panel(tiles, order),
        population_strip(order, colours),
        dendrogram_panel(clustering, colours, labels),
        design = "A##\nBCD", # The track spans only the heatmap column it shares an axis with
        widths = c(1L, 0.009, 0.1),
        heights = c(1L, 2.25)
    )

    ggsave(
        filename = paste0(populations[["A"]], "_", populations[["B"]], ".dcms.main.png"),
        plot = figure,
        dpi = 600L,
        height = figure_height(length(stats)),
        width = 6.75,
        bg = "white"
    )
    return()
}

walk(levels(dcms$COMPARISON), draw_figure)
