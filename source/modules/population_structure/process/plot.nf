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
