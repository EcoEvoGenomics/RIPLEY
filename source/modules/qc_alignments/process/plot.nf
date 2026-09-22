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

