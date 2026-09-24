process ADMIXTURE {

    label "ADMIXTURE"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    each(K)

    output:
    tuple path("${bed.simpleName}.k${K}.out"), path("${bed.simpleName}.k${K}.P"), path("${bed.simpleName}.k${K}.Q"), emit: data
    path("${bed.simpleName}.k${K}.alleles"), emit: alleles
    path("${bed.simpleName}.k${K}.clust"), emit: clust
    tuple val(K), env("cv_error"), emit: error

    script:
    """
    awk '{\$1="0";print \$0}' ${bed.simpleName}.bim > ${bed.simpleName}.bim.tmp
    mv ${bed.simpleName}.bim.tmp ${bed.simpleName}.bim

    admixture --cv -j${task.cpus} ${bed.simpleName}.bed ${K} > ${bed.simpleName}.k${K}.out
    mv ${bed.simpleName}.${K}.P ${bed.simpleName}.k${K}.P
    mv ${bed.simpleName}.${K}.Q ${bed.simpleName}.k${K}.Q

    awk 'NR==FNR {snp[FNR]=\$2; a1[FNR]=\$5; a2[FNR]=\$6; next} {print snp[FNR], a1[FNR], a2[FNR], \$0}' ${bim} ${bed.simpleName}.k${K}.P > ${bed.simpleName}.k${K}.alleles
    awk 'NR==FNR {id[FNR]=\$1; next} {print ${K}, id[FNR], \$0}' ${fam} ${bed.simpleName}.k${K}.Q > ${bed.simpleName}.k${K}.clust

    cv_error=\$(grep "CV error (K=" "${bed.simpleName}.k${K}.out" | awk '{print \$NF}')
    cv_error=\${cv_error:-NA}
    export cv_error
    """
}

process ADMIXTURE_PARENTAL_CENSUSES {

    // Expected format of input clustfile:
    // K  SAMPLE_ID  q(K1)  q(K2)  ...  q(K)

    label "RDATA"

    input:
    path(clustfile)
    val(parental_threshold)

    output:
    path("p*.list")

    script:
    """
    #!/usr/bin/env Rscript

    threshold <- as.numeric("${parental_threshold}")

    clust <- data.table::fread("${clustfile}", header = FALSE)
    ids <- as.character(clust[[2]])
    qmatrix <- as.matrix(clust[, 3:ncol(clust)])
    k <- ncol(qmatrix)

    censuses <- lapply(seq_len(k), function(i) ids[which(qmatrix[, i] > threshold)])
    empty <- which(lengths(censuses) == 0)

    if (length(empty) > 0) {
        stop(
            "No sample is assigned to cluster(s) ",
            paste(paste("p", empty, sep = ""), collapse = ", ")
        )
    }

    for (i in seq_len(k)) {
        writeLines(censuses[[i]], paste("p", i, ".list", sep = ""))
    }
    """
}

process IDENTIFY_AIMS {

    // Using between-population allele frequency variances to find AIMs,
    // inspired by https://doi.org/10.3389/fgene.2019.00043

    // Expected format of each input frqfile (vcftools --freq, header skipped):
    // CHROM  POS  N_ALLELES  N_CHR  REF:p(REF)  ALT:p(ALT)

    label "RDATA"

    input:
    path(frqfiles)
    path(sitekeys)
    val(variance_threshold)
    val(k)

    output:
    path("*.alleles"), optional: true, emit: alleles
    path("*.snpids"), optional: true, emit: snpids

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    threshold <- as.numeric("${variance_threshold}")
    k <- as.integer("${k}")

    sitekeys <- fread("${sitekeys}", header = FALSE, col.names = c("CHROM", "POS", "LOC"))

    # Every population is a subset of one VCF, so REF is shared and its frequency is comparable
    read_frequencies <- function(path, population) {
        frq <- fread(path, skip = 1, header = FALSE, fill = TRUE, sep = "\\t")
        setnames(frq, 1:6, c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF", "ALT"))
        frq <- frq[N_ALLELES == 2 & N_CHR > 0]
        out <- frq[, list(
            CHROM = CHROM,
            POS = POS,
            A1 = sub(":.*", "", REF),
            A2 = sub(":.*", "", ALT),
            P = as.numeric(sub(".*:", "", REF))
        )]
        setnames(out, "P", paste("P", population, sep = ""))
        out
    }

    frqfiles <- list.files(pattern = "[.]frq\$")
    populations <- as.integer(sub(".*_p([0-9]+)[.]frq\$", "\\\\1", frqfiles))
    frqfiles <- frqfiles[order(populations)]
    populations <- sort(populations)

    frequencies <- Reduce(
        function(x, y) merge(x, y, by = c("CHROM", "POS", "A1", "A2")),
        Map(read_frequencies, frqfiles, populations)
    )
    frequencies <- merge(frequencies, sitekeys, by = c("CHROM", "POS"))

    for (a in seq_along(populations)) {
        if (a == length(populations)) break
        for (b in seq(a + 1, length(populations))) {
            pop_a <- populations[a]
            pop_b <- populations[b]
            p_a <- frequencies[[paste("P", pop_a, sep = "")]]
            p_b <- frequencies[[paste("P", pop_b, sep = "")]]

            # Variance of a pair of frequencies, so caps at 0.5
            variances <- apply(cbind(p_a, p_b), 1, var)
            aims <- frequencies[variances >= threshold]
            if (nrow(aims) == 0) next

            q_a <- aims[[paste("P", pop_a, sep = "")]]
            q_b <- aims[[paste("P", pop_b, sep = "")]]
            aimtable <- data.table(
                LOC = aims[["LOC"]],
                P1 = paste("p", pop_a, sep = ""),
                A1 = fifelse(q_a > q_b, aims[["A1"]], aims[["A2"]]),
                P1_FREQA1 = fifelse(q_a > q_b, q_a, 1 - q_a),
                P2 = paste("p", pop_b, sep = ""),
                A2 = fifelse(q_b > q_a, aims[["A1"]], aims[["A2"]]),
                P2_FREQA2 = fifelse(q_b > q_a, q_b, 1 - q_b)
            )

            stem <- paste("aims_k", k, "_p", pop_a, "p", pop_b, sep = "")
            write.table(
                aimtable,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                file = paste(stem, ".alleles", sep = "")
            )
            writeLines(aimtable[["LOC"]], paste(stem, ".snpids", sep = ""))
        }
    }
    """
}

process CALCULATE_AIM_HIHET {

    label "RDATA"

    input:
    tuple path(allelefile), path(genotable)

    output:
    path("${genotable.simpleName}.aims"), emit: aims
    path("${genotable.simpleName}.hihet"), emit: hihet

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    al <- read.table("${allelefile.toString()}", col.names = c("LOC", "P1", "A1", "P1_FREQA1", "P2", "A2", "P2_FREQA2"))
    gt <- read.table("${genotable.toString()}", header = TRUE) |> rename(LOC = ID)

    data <- left_join(al, gt, by = "LOC") |> mutate(across(c(P1_FREQA1, P2_FREQA2), ~ sprintf("%.3f", .)))
    write.table(data, file = "${genotable.simpleName}.aims", quote = FALSE, row.names = FALSE, sep = "\\t")

    data <- data |>
        select(-c("P1_FREQA1", "P2_FREQA2")) |>
        mutate(COMPARISON = paste(P1, P2, sep = "-")) |>
        select(-c("P1", "P2")) |>
        pivot_longer(
            -c("LOC", "COMPARISON", "A1", "A2"),
            names_to = "ID",
            values_to = "GT"
        ) |>
        separate(GT, into = c("GT1", "GT2"), sep = "[/|]") |>
        mutate(
            called = !(GT1 == "." | GT2 == "."),
            nA1 = (GT1 == A1) + (GT2 == A1),
            nA2 = (GT1 == A2) + (GT2 == A2),
            het = (GT1 == A1 & GT2 == A2) | (GT1 == A2 & GT2 == A1)
        ) |>
        group_by(COMPARISON, ID) |>
        summarise(
            n_called = sum(called),
            n_missing = sum(!called),
            n_total = n_called + n_missing,
            n_A1 = sum(nA1[called]),
            n_A2 = sum(nA2[called]),
            n_het = sum(het[called]),
            HI = n_A1 / (n_A1 + n_A2),
            HET = n_het / n_called,
            MISS = n_missing / n_total
        ) |>
        select(COMPARISON, ID, HI, HET, MISS) |>
        mutate(across(c(HI, HET, MISS), ~ sprintf("%.3f", .)))
    
    write.table(data, file = "${genotable.simpleName}.hihet", quote = FALSE, row.names = FALSE, sep = "\\t")
    """
}
