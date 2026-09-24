process VCFTOOLS_FILTER_VCF {

    label "BCFTOOLS_VCFTOOLS"

    input:
    path(vcf)
    path(filter_flags)

    output:
    path("${vcf.simpleName}.vcf.gz")

    script:
    """
    echo "--gzvcf ${vcf}" > vcftools.args
    echo "--remove-indels" >> vcftools.args
    echo "--recode-INFO-all" >> vcftools.args
    echo "--recode" >> vcftools.args
    echo "--stdout" >> vcftools.args
    cat ${filter_flags} >> vcftools.args

    cat vcftools.args \
        | xargs vcftools \
        | bcftools view \
            --threads ${task.cpus} \
            --exclude 'ALT="*" || TYPE!="snp"' \
            --output-type z --output ${vcf.simpleName}_tmp.vcf.gz
    mv ${vcf.simpleName}_tmp.vcf.gz ${vcf.simpleName}.vcf.gz
    """
}

