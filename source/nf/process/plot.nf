process PLOT_ADMIXTURE {

    label "RPLOT"

    input:
    path(rscript)
    path(admixture_clusts)
    val(k_min_error)
    path(sample_metadata)
    path(population_metadata)
    path(species_metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${admixture_clusts} ${k_min_error} ${sample_metadata} ${population_metadata} ${species_metadata}
    """
}

process PLOT_HIHET {

    label "RPLOT"

    input:
    path(rscript)
    tuple path(hihet), path(sample_metadata), path(population_metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${hihet} ${sample_metadata} ${population_metadata}
    """
}

process PLOT_PLINK_LD_DECAY {

    label "RPLOT"

    input:
    path(rscript)
    path(ld_decay)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${ld_decay}
    """
}

process PLOT_PLINK_PCA {

    label "RPLOT"

    input:
    path(rscript)
    path(eigenval)
    path(eigenvec)
    path(sample_metadata)
    path(population_metadata)
    path(species_metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${eigenval} ${eigenvec} ${sample_metadata} ${population_metadata} ${species_metadata}
    """
}

process PLOT_SAMTOOLS_CRAM_STATS {

    label "RPLOT"

    input:
    path(rscript)
    path(stats)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${stats}
    """
}

process PLOT_SAMTOOLS_CRAM_STATS_POPWISE {

    label "RPLOT"

    input:
    path(rscript)
    path(stats)
    path(population_metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${stats} ${population_metadata}
    """
}

process PLOT_SAMTOOLS_BEDCOV {

    // bedcov_reference fixes one depth scale across all subsets so popwise plots are comparable

    label "RPLOT"

    input:
    path(rscript)
    each(bedcov)
    val(chrom_string)
    path(chrom_labels)
    path(bedcov_reference, stageAs: "reference.bedcov")

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${bedcov} ${chrom_string} ${chrom_labels} ${bedcov_reference}
    """
}

process PLOT_VCFTOOLS_RELATEDNESS {

    label "RPLOT"

    input:
    path(rscript)
    path(relatedness)
    path(sample_metadata)
    path(population_metadata)
    path(species_metadata)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${relatedness} ${sample_metadata} ${population_metadata} ${species_metadata}
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RPLOT"

    input:
    path(rscript)
    each(snpden)
    val(chrom_string)
    path(chrom_labels)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${snpden} ${chrom_string} ${chrom_labels}
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RPLOT"

    input:
    path(rscript)
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
    Rscript ${rscript} ${frq} ${idepth} ${imiss} ${ldepth_mean} ${lqual} ${lmiss} ${het} ${hwe} 
    """
}

process PLOT_VCFTOOLS_VCF_STATS_POPWISE {

    label "RPLOT"

    input:
    path(rscript)
    path(stats, stageAs: "stats/*")
    path(population_metadata)
    
    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} stats/ ${population_metadata}
    """
}

process PLOT_VCFTOOLS_PAIRWISE_MEAN_FST {

    label "RPLOT"

    input:
    path(rscript)
    path(means)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} ${means}
    """
}
