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
