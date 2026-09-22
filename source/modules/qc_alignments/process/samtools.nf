process SAMTOOLS_STATS {

    // Keeps only the SN summary section as key-value pairs

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai)

    output:
    path("${cram.simpleName}.stats")

    script:
    """
    samtools stats ${cram} \
        | awk -F '\\t' 'BEGIN { OFS = "\\t" } \$1 == "SN" { sub(/:\$/, "", \$2); print \$2, \$3 }' \
        > ${cram.simpleName}.stats
    """
}

process SAMTOOLS_BEDCOV {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai), path(windows_bed)

    output:
    path("${cram.simpleName}.bedcov")

    script:
    """
    samtools bedcov ${windows_bed} ${cram} > ${cram.simpleName}.bedcov
    """
}
