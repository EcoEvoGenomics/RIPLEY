process PLOT_DCMS {

    label "RPLOT"

    input:
    path(rscript)
    path(dcms_tables, stageAs: "dcms/*")
    path(population_metadata)
    path(chrom_labels)

    output:
    path("*.png")

    script:
    """
    Rscript ${rscript} dcms/ ${population_metadata} ${chrom_labels}
    """
}
