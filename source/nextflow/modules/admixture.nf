process ADMIXTURE {

    label "ADMIXTURE"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    each(K)

    output:
    tuple path("${bed.simpleName}.k${K}.out"), path("${bed.simpleName}.k${K}.Q"), path("${bed.simpleName}.k${K}.P"), emit: data
    path("${bed.simpleName}.k${K}.P"), emit: pfile
    path("${bed.simpleName}.clust"), emit: clust
    tuple val(K), env("cv_error"), emit: error

    script:
    """
    awk '{\$1="0";print \$0}' ${bed.simpleName}.bim > ${bed.simpleName}.bim.tmp
    mv ${bed.simpleName}.bim.tmp ${bed.simpleName}.bim

    admixture --cv -j${task.cpus} ${bed.simpleName}.bed ${K} > ${bed.simpleName}.k${K}.out
    mv ${bed.simpleName}.${K}.P ${bed.simpleName}.k${K}.P
    mv ${bed.simpleName}.${K}.Q ${bed.simpleName}.k${K}.Q

    awk 'NR==FNR {f2[FNR]=\$1; next} {print ${K}, f2[FNR], \$0}' ${fam} ${bed.simpleName}.k${K}.Q > ${bed.simpleName}.clust

    cv_error=\$(grep "CV error (K=" "${bed.simpleName}.k${K}.out" | awk '{print \$NF}')
    cv_error=\${cv_error:-NA}
    export cv_error
    """
}

process ADMIXTURE_AIMS {

    // Using .P-file between-column variances to find AIMs,
    // inspired by https://doi.org/10.3389/fgene.2019.00043
    
    label "RBASE"

    input:
    each(pfile)
    path(snplist)
    val(variance_threshold)

    output:
    path("*.list"), optional: true

    script:
    """
    #!/usr/bin/env Rscript
    snp_ids <- read.table("${snplist.toString()}")\$V1
    p_table <- read.table("${pfile.toString()}")
    k <- ncol(p_table)
    for (i in seq_len(k)) {
        if (i == k) break
        for (j in seq(i + 1, k)) {
            ij_vars <- apply(p_table[c(i, j)], 1, \\(x) var(x))
            ij_aims <- snp_ids[which(ij_vars > ${variance_threshold})]
            writeLines(ij_aims, paste("aims_k", k, "_p", i, "p", j, ".list", sep = ""))
        }
    }
    """

}
