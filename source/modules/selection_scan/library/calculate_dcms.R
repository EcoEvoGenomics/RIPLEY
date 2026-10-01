# De-correlated composite of multiple signals (DCMS)
# Ma et al. 2015 https://doi.org/10.1038/hdy.2015.42
#
# This implementation after a script by Yuheng Sun.

options(scipen = 999L)

# Enforced upstream in recast_for_dcms.R
dcms_layout <- c(
    "CHR",
    "WINDOW_START",
    "WINDOW_END",
    "POP_A",
    "POP_B",
    "STAT_VALUE",
    "DIRECTION"
)

# Enforced upstream in recast_for_dcms.R
valid_directions <- c(
    "low",
    "high",
    "two_tailed"
)

read_stat_table <- function(path) {
    stat <- sub("\\.dcms\\.tsv$", "", basename(path))
    stat <- tail(strsplit(stat, ".", fixed = TRUE)[[1L]], 1L)
    table <- read.delim(
        path,
        header = TRUE,
        na.strings = c("NA", "nan", ""),
        colClasses = c(CHR = "character", POP_A = "character", POP_B = "character")
    )

    layout_valid <- identical(colnames(table), dcms_layout)
    if (! layout_valid) {
        stop(basename(path), "has invalid layout.", call. = FALSE)
    }

    direction <- unique(table$DIRECTION)
    direction_valid <- direction %in% valid_directions
    direction_unambiguous <- length(direction) == 1L
    if (! direction_valid || ! direction_unambiguous) {
        stop(basename(path), " must carry one valid DIRECTION.", call. = FALSE)
    }

    table$STAT <- stat
    return(table)
}

collate_stat_tables_for_pair <- function(all_tables, pop_a, pop_b) {

    shared_table_keys <- c("CHR", "WINDOW_START", "WINDOW_END")
    collated_columns <- list()
    directions <- character(0L)

    for (table in all_tables) {

        stat <- table$STAT[1L]
        stat_is_pairwise <- any(!is.na(table$POP_B))

        if (! stat_is_pairwise) {
            row_pop_is_a <- table$POP_A == pop_a
            row_pop_is_b <- table$POP_A == pop_b
            rows_with_a <- table[row_pop_is_a, ]
            rows_with_b <- table[row_pop_is_b, ]
            to_collate <- list(rows_with_a, rows_with_b)
            names(to_collate) <- c(paste(stat, "A", sep = "_"), paste(stat, "B", sep = "_"))
        } else {
            row_pop_is_both <- table$POP_A == pop_a & table$POP_B == pop_b
            row_pop_is_htob <- table$POP_A == pop_b & table$POP_B == pop_a
            rows_with_both <- table[row_pop_is_both | row_pop_is_htob, ]
            to_collate <- list(rows_with_both)
            names(to_collate) <- stat
        }

        # Accumulate stats by row to collect chromosomes from multiple files before joining columns
        for (column_name in names(to_collate)) {
            rows <- to_collate[[column_name]]
            if (nrow(rows) == 0L) {
                next
            }
            directions[column_name] <- rows$DIRECTION[1L]
            rows <- rows[, c(shared_table_keys, "STAT_VALUE")]
            colnames(rows)[colnames(rows) == "STAT_VALUE"] <- column_name
            collated_columns[[column_name]] <- rbind(collated_columns[[column_name]], rows)
        }

    }

    # Windows scans disagree on are not joined (e.g. rehh windows overhang contig end, pixy not)
    joined <- NULL
    for (column_name in names(collated_columns)) {
        joined <- if (is.null(joined)) {
            collated_columns[[column_name]]
        } else {
            merge(joined, collated_columns[[column_name]], by = shared_table_keys, all = FALSE)
        }
    }

    return(list(table = joined, directions = directions))
}

apply_value_directionality <- function(value, direction) {
    directional <- switch(direction,
        high = value,
        low = -value,
        two_tailed = abs(value),
        stop("Invalid DIRECTION '", direction, "'.", call. = FALSE)
    )
    return(directional)
}

calculate_fractional_rank_pval <- function(value) {
    observed <- !is.na(value)
    p <- rep(NA_real_, length(value))
    if (sum(observed) < 2L) {
        return(p)
    }
    fractional_rank <- rank(value[observed]) / (sum(observed) + 1.0)
    normalised <- scale(qnorm(fractional_rank))
    p[observed] <- pnorm(normalised, lower.tail = FALSE)
    return(p)
}

calculate_dcms <- function(collated_stats, directions) {
    stopifnot(all(colnames(collated_stats) %in% names(directions)))

    direction_map <- directions[colnames(collated_stats)]
    directed <- Map(apply_value_directionality, collated_stats, direction_map)
    directed <- as.data.frame(directed)

    correlations <- cor(directed, method = "spearman", use = "pairwise.complete.obs")
    correlations[is.na(correlations)] <- 0L
    weights <- colSums(abs(correlations))

    p_values <- as.data.frame(lapply(directed, calculate_fractional_rank_pval))
    colnames(p_values) <- colnames(collated_stats)
    logodds_scores <- log((1.0 - p_values) / p_values)
    dcms_values <- rowSums(sweep(as.matrix(logodds_scores), 2L, weights, "/"))

    # An NA value for any stat deliberately yields NA for whole window: rowSums(na.rm = FALSE)
    return(list(p = p_values, dcms = dcms_values))
}

main <- function() {

    table_paths <- sort(list.files("dcms", pattern = "\\.dcms\\.tsv$", full.names = TRUE))
    if (length(table_paths) == 0L) {
        stop(
            "No files were staged as ./dcms/<stat>.dcms.tsv.",
            call. = FALSE
        )
    }

    tables <- lapply(table_paths, read_stat_table)

    pairwise_rows <- do.call(rbind, lapply(tables, \(x) x[!is.na(x$POP_B), c("POP_A", "POP_B")]))
    if (is.null(pairwise_rows) || nrow(pairwise_rows) == 0L) {
        stop(
            "At least one pairwise statistic required to define population pairs for DCMS.",
            call. = FALSE
        )
    }

    population_pairs <- unique(t(apply(pairwise_rows, 1L, sort)))

    for (pair_index in seq_len(nrow(population_pairs))) {

        pop_a <- population_pairs[pair_index, 1L]
        pop_b <- population_pairs[pair_index, 2L]

        collated <- collate_stat_tables_for_pair(tables, pop_a, pop_b)
        keys <- c("CHR", "WINDOW_START", "WINDOW_END")

        value_columns <- setdiff(colnames(collated$table), keys)
        if (length(value_columns) < 2L) {
            stop(
                "Pair ", pop_a, "-", pop_b, " retained only ", length(value_columns),
                " statistic(s); DCMS requires at least two.",
                call. = FALSE
            )
        }

        # Only windows where all statistics are measured are comparable; drop incomplete
        complete <- stats::complete.cases(collated$table[, value_columns, drop = FALSE])
        if (any(!complete)) {
            collated$table <- collated$table[complete, , drop = FALSE]
        }

        if (nrow(collated$table) < 2L) {
            stop(
                "Pair ", pop_a, "-", pop_b, " retained ", nrow(collated$table),
                " window(s) measured by every statistic; DCMS ranks across at least two.",
                call. = FALSE
            )
        }

        dcms <- calculate_dcms(
            collated$table[, value_columns, drop = FALSE],
            collated$directions
        )

        # Values are published as in the original tables, directionality is only applied internally
        output <- collated$table[, keys]
        output[value_columns] <- collated$table[, value_columns, drop = FALSE]
        output[paste0(value_columns, "_P")] <- dcms$p
        output$DCMS <- dcms$dcms
        output <- output[order(output$CHR, output$WINDOW_START), ]

        write.table(output, file = paste0(pop_a, "_", pop_b, ".dcms.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    }

    return()

}

main()
