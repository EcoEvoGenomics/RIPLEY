process PIXY_STATS {

    label "PIXY"

    input:
    tuple path(vcf), path(csi), path(windows_bed)
    path(popmap)
    each(stat)

    output:
    path("${vcf.simpleName}.${stat}.tsv")

    script:
    """
    # windows_bed is genome-wide, but each task should see only one chromosome
    awk -v OFS="\\t" '\$1 == "${vcf.simpleName}"' ${windows_bed} > chrom_windows.bed
    if [ ! -s chrom_windows.bed ]; then
        echo "No windows on chromosome ${vcf.simpleName} in ${windows_bed}." >&2
        exit 1
    fi

    pixy --stats ${stat} \
        --vcf ${vcf} \
        --populations ${popmap} \
        --bed_file chrom_windows.bed \
        --n_cores ${task.cpus} \
        --output_prefix pixy
    mv pixy_${stat}.txt ${vcf.simpleName}.${stat}.tsv
    """
}

process RECAST_PIXY_FOR_DCMS {

    label "RBASE"

    input:
    tuple path(tsv), val(stat)
    path(recast_rscript)

    output:
    path("${tsv.simpleName}.${stat}.dcms.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    options(scipen = 999)
    source("${recast_rscript.toString()}")

    write_dcms_formatted(
        recast_pixy_for_dcms("${tsv.toString()}", "${stat}"),
        "${tsv.simpleName}.${stat}.dcms.tsv"
    )
    """
}
