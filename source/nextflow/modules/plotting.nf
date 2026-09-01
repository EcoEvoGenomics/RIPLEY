process PLOT_PLINK_LD_DECAY {

    label "RPLOT"

    input:
    path(plot_script)
    path(ld_decay)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${ld_decay}
    """
}

process PLOT_PLINK_PCA {

    label "RPLOT"

    input:
    path(plot_script)
    path(eigenval)
    path(eigenvec)
    path(metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${eigenval} ${eigenvec} ${metadata}
    """
}

process PLOT_REHH_XPEHH {

    label "RPLOT"

    input:
    path(plot_script)
    path(scans)
    path(cands)
    path(gff)
    path(chrom_labels)
    val(cand_pval)

    output:
    path("*.png"), emit: mainplot
    path("candidate_regions.csv"), emit: candregions
    path("**/*.png"), emit: candplots
    path("**/*.gff"), emit: candgenes

    script:
    """
    Rscript ${plot_script} ${scans} ${cands} ${gff} ${chrom_labels} ${cand_pval} 170 170 20
    """
}

process PLOT_VCFTOOLS_RELATEDNESS {

    label "RPLOT"

    input:
    path(plot_script)
    path(relatedness)
    path(metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${relatedness} ${metadata}
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RPLOT"

    input:
    path(plot_script)
    each(snpden)
    val(chrom_string)
    path(chrom_labels)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${snpden} ${chrom_string} ${chrom_labels}
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RBASE"

    input:
    path(plot_script)
    tuple \
        path(frq),
        path(idepth),
        path(imiss),
        path(ldepth_mean),
        path(lqual),
        path(lmiss),
        path(het),
        path(hwe)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${frq} ${idepth} ${imiss} ${ldepth_mean} ${lqual} ${lmiss} ${het} ${hwe} 
    """
}

process PLOT_VCFTOOLS_PAIRWISE_MEAN_FST {

    label "RPLOT"

    input:
    path(plot_script)
    path(means)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${means}
    """
}
