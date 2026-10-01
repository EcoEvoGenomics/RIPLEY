# Tables for calculate_dcms.R must pass through here.

#---- COMMON -------------------------------------------------------------------
# Checked downstream in calculate_dcms.R
dcms_layout <- c(
    "CHR",
    "WINDOW_START",
    "WINDOW_END",
    "POP_A",
    "POP_B",
    "STAT_VALUE",
    "DIRECTION"
)

dcms_formatted <- function(chr, start, end, pop_a, pop_b, value, direction) {
    table <- data.frame(
        CHR = chr,
        WINDOW_START = start,
        WINDOW_END = end,
        POP_A = pop_a,
        POP_B = pop_b,
        STAT_VALUE = value,
        DIRECTION = direction,
        stringsAsFactors = FALSE
    )
    stopifnot(identical(colnames(table), dcms_layout))
    return(table)
}

write_dcms_formatted <- function(table, file) {
    write.table(
        table,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    return()
}

#---- PIXY ---------------------------------------------------------------------
pixy_schemas <- list(
    pi = list(
        value = "avg_pi", pops = c("pop"), direction = "low"
    ),
    watterson_theta = list(
        value = "avg_watterson_theta", pops = c("pop"), direction = "low"
    ),
    tajima_d = list(
        value = "tajima_d", pops = c("pop"), direction = "two_tailed"
    ),
    dxy = list(
        value = "avg_dxy", pops = c("pop1", "pop2"), direction = "high"
    ),
    fst = list(
        value = "avg_wc_fst", pops = c("pop1", "pop2"), direction = "high"
    )
)

recast_pixy_for_dcms <- function(path, stat) {

    stat_valid <- stat %in% names(pixy_schemas)
    if (! stat_valid) {
        stop("No DCMS schema for pixy statistic '", stat, "'.", call. = FALSE)
    }
    schema <- pixy_schemas[[stat]]

    pixy <- read.delim(path, header = TRUE, na.strings = c("NA", "nan", ""))

    required_cols <- c(
        "chromosome",
        "window_pos_1",
        "window_pos_2",
        schema$pops,
        schema$value
    )
    missing_cols <- setdiff(required_cols, colnames(pixy))

    if (length(missing_cols) > 0L) {
        stop(basename(path), " lacks expected column(s) for '", stat, "': ",
             paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }

    stat_is_pairwise <- length(schema$pops) > 1L
    if (stat_is_pairwise) {
        pop_b <- pixy[[schema$pops[2L]]]
    } else {
        pop_b <- NA_character_
    }

    formatted <- dcms_formatted(
        chr = pixy$chromosome,
        start = pixy$window_pos_1,
        end = pixy$window_pos_2,
        pop_a = pixy[[schema$pops[1L]]],
        pop_b = pop_b,
        value = pixy[[schema$value]],
        direction = schema$direction
    )
    return(formatted)
}

#---- REHH ---------------------------------------------------------------------
# 'threshold' is arbitrary because p-values are calculated by DCMS downstream,
# but calc_candidate_regions rejects zero; ignore_sign so that the population on
# the negative side of a comparison is not zeroed away before window averaging.

bin_rehh_output <- function(scan, window_size, step_size, min_sites) {

    # rehh's 'overlap' is a stride that cannot equal the window size.
    stride <- if (step_size >= window_size) 0L else as.integer(step_size)

    binned <- rehh::calc_candidate_regions(
        scan = scan,
        threshold = 1e-9,
        ignore_sign = TRUE,
        window_size = window_size,
        overlap = stride,
        min_n_mrk = min_sites,
        min_n_extr_mrk = 0L,
        join_neighbors = FALSE
    )
    return(binned)
}

# rehh windows are half-open from 0, but shared DCMS layout is 1-based inclusive
recast_rehh_bins_for_dcms <- function(windows, pop_a, pop_b) {
    formatted <- dcms_formatted(
        chr = windows$CHR,
        start = windows$START + 1L,
        end = windows$END,
        pop_a = pop_a,
        pop_b = pop_b,
        value = windows$MEAN_MRK,
        direction = "high"
    )
    return(formatted)
}
