process PIXY_STATS {

    label "PIXY"

    input:
    tuple path(vcf), path(csi)
    path(popmap)
    each(stat)
    val(window_size)

    output:
    path("${vcf.simpleName}.${stat}.tsv")

    script:
    """
    pixy --stats ${stat} \
        --vcf ${vcf} \
        --populations ${popmap} \
        --window_size ${window_size} \
        --n_cores ${task.cpus} \
        --output_prefix pixy
    mv pixy_${stat}.txt ${vcf.simpleName}.${stat}.tsv
    """
}
