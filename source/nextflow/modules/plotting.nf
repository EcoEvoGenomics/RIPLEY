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

process PLOT_REHH_XPEHH {

    label "RPLOT"

    input:
    path(plot_script)
    path(scans)
    path(cands)
    path(gff)
    path(chrom_conversions)
    val(cand_pval)

    output:
    path("*.png"), emit: mainplot
    path("candidate_regions.csv"), emit: candregions
    path("**/*.png"), emit: candplots
    path("**/*.gff"), emit: candgenes

    script:
    """
    Rscript ${plot_script} ${scans} ${cands} ${gff} ${chrom_conversions} ${cand_pval} 170 170 20
    """
}

process PLOT_VCFTOOLS_RELATEDNESS {

    label "RPLOT"

    input:
    path(plot_script)
    path(relatedness)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${relatedness}
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RPLOT"

    input:
    path(plot_script)
    path(snpden)
    val(chrom_string)
    path(chrom_conversions)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${snpden} ${chrom_string} ${chrom_conversions}
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RBASE"

    input:
    path(plot_script)
    path(frq)
    path(idepth)
    path(imiss)
    path(ldepth_mean)
    path(lqual)
    path(lmiss)
    path(het)
    path(hwe)

    output:
    path("*.png")

    script:
    """
    Rscript ${plot_script} ${frq} ${idepth} ${imiss} ${ldepth_mean} ${lqual} ${lmiss} ${het} ${hwe} 
    """
}
