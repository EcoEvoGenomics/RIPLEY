process PLOT_PLINK_LD_DECAY {

    label "RPLOT"

    input:
    path(ld_decay_bins)
    val(ld_window_kb)

    output:
    path("${ld_decay_bins}.png")

    script:
    """
    #!/usr/bin/env Rscript
    library(ggplot2)
    tbl <- read.table("${ld_decay_bins.toString()}", header = TRUE)
    ordered_chrs <- unique(tbl\$CHR)
    tbl\$CHR <- factor(tbl\$CHR, levels = ordered_chrs)
    chr_count <- length(unique(tbl\$CHR))
    png("${ld_decay_bins}.png", height = chr_count / 2, width = chr_count / 2, units = "in", res = 600)
    tbl |>
        ggplot(aes(x = DIST, y = AVG_R2)) +
        geom_point(alpha = 0.1, stroke = NA) +
        facet_wrap(vars(CHR)) +
        scale_y_continuous(name = bquote("Linkage diseqiulibrium, mean" ~ R^2), n.breaks = 3, limits = c(0, 1), expand = FALSE) +
        scale_x_continuous(name = "Distance (bp)", n.breaks = 3) +
        theme_bw() +
        theme(
            panel.spacing = unit(10, "mm")
        )
    dev.off()
    """
}

process PLOT_REHH_XPEHH {

    label "RPLOT"

    input:
    path(parsescript)
    path(scans)
    path(cands)
    path(gff)
    path(chrom_renames)
    val(cand_pval)

    output:
    path("*.png"), emit: mainplot
    path("candidate_regions.csv"), emit: candregions
    path("**/*.png"), emit: candplots
    path("**/*.gff"), emit: candgenes

    script:
    """
    Rscript ${parsescript} ${scans} ${cands} ${gff} ${chrom_renames} ${cand_pval} 170 170 20
    """
}

process PLOT_VCFTOOLS_RELATEDNESS {

    label "RPLOT"

    input:
    path(relatedness)

    output:
    path("${relatedness.simpleName}.png")

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(pheatmap)

    tbl <- read.table("${relatedness.toString()}", header = TRUE) |>
        select(INDV1, INDV2, RELATEDNESS_PHI) |>
        pivot_wider(names_from = INDV2, values_from = RELATEDNESS_PHI) |>
        column_to_rownames("INDV1")

    second_degree_relatedness <- 1/(2^(7/2))
    highlights <- as.matrix(tbl)
    highlights <- formatC(highlights, format = "g", digits = 2)
    highlights[is.na(highlights) | highlights <= second_degree_relatedness] <- ""

    plt <- pheatmap(
        mat = tbl,
        color = viridisLite::inferno(100),
        clustering_method = "ward.D2",
        border_color = NA,
        legend = FALSE,
        display_numbers = highlights,
        number_color = "red",
        fontsize_number = 5
    )

    n_inds <- nrow(tbl)
    inches_per_ind <- 0.25

    ggsave(
        plot = plt,
        filename = "${relatedness.simpleName}.png",
        dpi = 600,
        height = n_inds * inches_per_ind,
        width = n_inds * inches_per_ind,
        bg = "white"
    )
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RPLOT"

    input:
    path(snpden)
    path(chrom_conversions)
    val(plot_chroms)

    output:
    path("${snpden}.png")

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    chroms <- strsplit("${plot_chroms}", ",")[[1]]

    tbl <- read.table("${snpden}", header = TRUE)
    tbl <- filter(tbl, CHROM %in% chroms)

    chr_conversions <- read.table("${chrom_conversions}")
    chr_conversions <- filter(chr_conversions, V1 %in% chroms)
    renamed_chrs <- chr_conversions\$V1
    names(renamed_chrs) <- chr_conversions\$V2

    bin_size <- tbl\$BIN_START[2] - tbl\$BIN_START[1]
    xmax <- max(tbl\$BIN_START) + bin_size

    plt <- tbl |>
    ggplot(
        aes(
        x = BIN_START,
        y = factor(CHROM, levels = rev(chroms)),
        fill = log10(SNP_COUNT)
        )
    ) +
    geom_tile(height = 0.75) +
    scale_fill_viridis_c(
        name = expression(log[10] ~ (SNPs)),
        option = "magma"
    ) +
    scale_x_continuous(
        guide = guide_axis(cap = TRUE),
        expand = expansion(add = 0),
        limits = c((xmax * -0.015), (xmax * 1.03)),
        breaks = seq(from = 0, to = xmax, length.out = 3),
        labels = scales::label_number(
        accuracy = 1,
        scale    = 1 / 1e6,
        suffix   = " mbp"
        )
    ) +
    scale_y_discrete(
        expand = c(0.025),
        labels = rev(names(renamed_chrs))
    ) +
    theme_void() +
    theme(
        axis.line.x = element_line(colour = "black", linewidth = 0.1),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, hjust = 1),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.1),
        legend.text = element_text(size = 6),
        legend.ticks = element_blank(),
        legend.title = element_text(size = 6, face = "bold", hjust = 1),
        legend.key.width  = unit(0.3, "cm"),
        legend.key.height = unit(0.4, "cm")
    )

    n_chroms <- length(chroms)
    inches_per_chrom <- 0.15

    ggsave(
        plot = plt,
        filename = "${snpden}.png",
        dpi = 600,
        height = n_chroms * inches_per_chrom,
        width = 6.75,
        bg = "white"
    )
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RBASE"

    input:
    path(vcf)
    path(frq)
    path(idepth)
    path(imiss)
    path(ldepth_mean)
    path(lqual)
    path(lmiss)
    path(het)
    path(hwe)

    output:
    path("${vcf.simpleName}.stats.pdf")

    script:
    """
    #!/usr/bin/env Rscript

    # Count greatest number of alleles in VCF to format frq table for MAF
    frq <- "${frq.toString()}"
    skip_header_line <- 1L

    field_counts <- count.fields(frq)
    field_counts <- field_counts[-seq_len(skip_header_line)]
    n_nonallele_cols <- 4
    n_alleles <- max(field_counts, na.rm = TRUE) - n_nonallele_cols
    col_names <- c("CHROM","POS","N_ALLELES","N_CHR", paste0("A", seq_len(n_alleles)))

    frq <- read.table(frq, skip = skip_header_line, header = FALSE, fill = TRUE, col.names = col_names, na.strings = c("", "NA"), stringsAsFactors = FALSE)
    
    frq\$MAF <- (frq[grep("^A", names(frq), value = TRUE)] |> apply(1, \\(x) min(x)))
    frq\$MAC <- round(frq\$MAF * frq\$N_CHR)
    idepth <- read.table("${idepth.toString()}", header = TRUE)
    imiss <- read.table("${imiss.toString()}", header = TRUE)
    ldepth_mean <- read.table("${ldepth_mean.toString()}", header = TRUE)
    lqual <- read.table("${lqual.toString()}", header = TRUE)
    lmiss <- read.table("${lmiss.toString()}", header = TRUE)
    het <- read.table("${het.toString()}", header = TRUE)
    hwe <- read.table("${hwe.toString()}", header = TRUE)

    pdf("${vcf.simpleName}.stats.pdf", width = 20, height = 20)
    par(mfrow = c(3, 3))
    hist(het\$F, main = "INDV INBREEDING COEFFICIENT (F)")
    hist(idepth\$MEAN_DEPTH, main = "INDV MEAN DEPTH")
    hist(imiss\$F_MISS, main = "INDV MISSINGNESS")
    hist(ldepth_mean\$MEAN_DEPTH[ldepth_mean\$MEAN_DEPTH <= quantile(ldepth_mean\$MEAN_DEPTH, 0.99)], main = "SITE MEAN DEPTH (Cutoff at 99th percentile)")
    hist(lqual\$QUAL[lqual\$QUAL <= quantile(lqual\$QUAL, 0.99)], main = "SITE QUALITY (Cutoff at 99th percentile)")
    hist(lmiss\$F_MISS, main = "SITE MISSINGNESS")
    hist(frq\$MAF, main = "SITE MINOR ALLELE FREQUENCY")
    hist(frq\$MAC, main = "SITE MINOR ALLELE COUNT", breaks = seq(min(frq\$MAC, na.rm = TRUE), max(frq\$MAC, na.rm = TRUE), length.out = max(frq\$MAC, na.rm = TRUE)))
    hist(-log10(hwe\$P_HWE), main = "SITE HARDY-WEINBERG TEST (-LOG10 P-VALUE)")
    dev.off()
    """
}
