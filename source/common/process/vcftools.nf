process VCFTOOLS_EXCLUDE_BED {

    label "VCFTOOLS"

    input:
    path(vcf)
    path(bed)

    output:
    path("${vcf.simpleName}.vcf.gz")

    script:
    """
    vcftools --gzvcf ${vcf} --exclude-bed ${bed} --recode --stdout \
        | gzip -c > ${vcf.simpleName}_tmp.vcf.gz
    mv ${vcf.simpleName}_tmp.vcf.gz ${vcf.simpleName}.vcf.gz
    """
}

