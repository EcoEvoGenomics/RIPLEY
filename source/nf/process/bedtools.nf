process BEDTOOLS_MAKEWINDOWS {

    label "BEDTOOLS"

    input:
    path(fai)
    val(size)
    val(step)

    output:
    path('windows.bed'), emit: bed_base_zero
    path('windows.txt'), emit: list_base_one

    script:
    """
    cut -f 1-2 ${fai} > contig_sizes.tsv
    bedtools makewindows -g contig_sizes.tsv -w ${size} -s ${step} > windows.bed
    awk '{print \$1":"\$2"-"\$3}' windows.bed | sed 's/:0-/:1-/g' > windows.txt
    """
}
