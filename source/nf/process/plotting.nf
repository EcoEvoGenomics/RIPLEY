process PLOT_PLINK_LD_DECAY {

    label "RPLOT"

    input:
    path(ld_decay)

    output:
    path("*.png")

    script:
    """
    Rscript ${projectDir}/../../R/plot_linkage_decay.R ${ld_decay}
    """
}

process PLOT_PLINK_PCA {

    label "RPLOT"

    input:
    path(eigenval)
    path(eigenvec)
    path(metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${projectDir}/../../R/plot_pca.R ${eigenval} ${eigenvec} ${metadata}
    """
}

process PLOT_REHH_XPEHH {

    label "RPLOT"

    input:
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
    Rscript ${projectDir}/../../R/plot_xpehh.R ${scans} ${cands} ${gff} ${chrom_labels} ${cand_pval} 170 170 20
    """
}

process PLOT_VCFTOOLS_RELATEDNESS {

    label "RPLOT"

    input:
    path(relatedness)
    path(metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${projectDir}/../../R/plot_kinship.R ${relatedness} ${metadata}
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RPLOT"

    input:
    each(snpden)
    val(chrom_string)
    path(chrom_labels)

    output:
    path("*.png")

    script:
    """
    Rscript ${projectDir}/../../R/plot_snpden.R ${snpden} ${chrom_string} ${chrom_labels}
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RPLOT"

    input:
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
    Rscript ${projectDir}/../../R/plot_vcfstats.R ${frq} ${idepth} ${imiss} ${ldepth_mean} ${lqual} ${lmiss} ${het} ${hwe} 
    """
}

process PLOT_VCFTOOLS_PAIRWISE_MEAN_FST {

    label "RPLOT"

    input:
    path(means)

    output:
    path("*.png")

    script:
    """
    Rscript ${projectDir}/../../R/plot_fst.R ${means}
    """
}
