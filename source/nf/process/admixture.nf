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

process ADMIXTURE_AIMS {

    // Using .P-file between-column variances to find AIMs,
    // inspired by https://doi.org/10.3389/fgene.2019.00043

    // Expected format of input allelefile:
    // SNP_ID  A1  A2  pA1(K1)   pA1(K2)   ...   pA1(K-1)   pA1(K)
    
    label "RDATA"

    input:
    path(allelefile)
    val(variance_threshold)

    output:
    path("*.alleles"), optional: true, emit: alleles
    path("*.snpids"), optional: true, emit: snpids

    script:
    """
    #!/usr/bin/env Rscript

    ptable <- read.table("${allelefile.toString()}")
    snpids <- ptable\$V1
    metacols <- 1:3
    nskip <- length(metacols)

    k <- ncol(ptable) - nskip

    for (i in seq_len(k)) {
        if (i == k) break
        for (j in seq(i + 1, k)) {
            p1 <- i + nskip
            p2 <- j + nskip
            vars <- apply(ptable[c(p1, p2)], 1, \\(x) var(x))
            aims <- snpids[which(vars >= ${variance_threshold})]
            aimtable <- ptable[which(ptable\$V1 %in% aims), c(metacols, p1, p2)]
            aimtable[6] <- ifelse(aimtable[[4]] > aimtable[[5]], aimtable[[2]], aimtable[[3]])
            aimtable[7] <- ifelse(aimtable[[4]] > aimtable[[5]], aimtable[[4]], 1 - aimtable[[4]])
            aimtable[8] <- ifelse(aimtable[[5]] > aimtable[[4]], aimtable[[2]], aimtable[[3]])
            aimtable[9] <- ifelse(aimtable[[5]] > aimtable[[4]], aimtable[[5]], 1 - aimtable[[5]])
            aimtable <- aimtable[, c(1, 6 : 9)]
            write.table(
                aimtable,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                file = paste("aims_k", k, "_p", i, "p", j, ".alleles", sep = "")
            )
            writeLines(aims, paste("aims_k", k, "_p", i, "p", j, ".snpids", sep = ""))
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

    al <- read.table("${allelefile.toString()}", col.names = c("LOC", "A1", "P1_FREQA1", "A2", "P2_FREQA2"))
    gt <- read.table("${genotable.toString()}", header = TRUE) |> rename(LOC = ID)

    data <- left_join(al, gt, by = "LOC") |> mutate(across(c(P1_FREQA1, P2_FREQA2), ~ sprintf("%.3f", .)))
    write.table(data, file = "${genotable.simpleName}.aims", quote = FALSE, row.names = FALSE, sep = "\\t")

    data <- data |>
        select(-c("P1_FREQA1", "P2_FREQA2")) |>
        pivot_longer(
            -c("LOC", "A1", "A2"),
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
        group_by(ID) |>
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
        select(ID, HI, HET, MISS) |>
        mutate(across(c(HI, HET, MISS), ~ sprintf("%.3f", .)))
    
    write.table(data, file = "${genotable.simpleName}.hihet", quote = FALSE, row.names = FALSE, sep = "\\t")
    """
}
