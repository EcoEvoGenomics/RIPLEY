process BEDTOOLS_MAKEWINDOWS {

    label "BEDTOOLS"

    input:
    path(fai)
    val(size)
    val(step)

    output:
    path('windows.bed'), emit: bed_base_zero
    path('windows.txt'), emit: regions_base_one

    script:
    """
    cut -f 1-2 ${fai} > contig_sizes.tsv
    bedtools makewindows -g contig_sizes.tsv -w ${size} -s ${step} > windows.bed
    awk -v OFS="\\t" '{print \$1"_"sprintf("%06d", ++n[\$1]), \$1":"\$2"-"\$3}' windows.bed | sed 's/:0-/:1-/g' > windows.txt
    """
}

process BEDTOOLS_INTERSECT_WINDOWS {

    // Drops variant-free windows from BEDTOOLS_MAKEWINDOWS

    label "BEDTOOLS"

    input:
    path(windows_bed)
    path(vcf)

    output:
    path("${vcf.simpleName}.windows.bed"), emit: bed_base_zero
    path("${vcf.simpleName}.windows.txt"), emit: regions_base_one

    script:
    """
    bedtools intersect -u -a ${windows_bed} -b ${vcf} > ${vcf.simpleName}.windows.bed
    awk -v OFS="\\t" '{print \$1"_"sprintf("%06d", ++n[\$1]), \$1":"\$2"-"\$3}' ${vcf.simpleName}.windows.bed \
        | sed 's/:0-/:1-/g' > ${vcf.simpleName}.windows.txt
    """
}
