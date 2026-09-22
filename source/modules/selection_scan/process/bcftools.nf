process BCFTOOLS_BCF_TO_VCF {

    label "BCFTOOLS"

    input:
    tuple path(bcf), path(csi)

    output:
    path("${bcf.simpleName}.vcf.gz")

    script:
    """
    bcftools view \
        --threads ${task.cpus} \
        --output-type z --output ${bcf.simpleName}.vcf.gz \
        ${bcf}
    """
}
