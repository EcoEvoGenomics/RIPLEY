process VCFTOOLS_FILTER_VARIANTS {

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

process VCFTOOLS_FILTER_VARIANTS_AND_INVARIANTS {

    label "BCFTOOLS_VCFTOOLS"

    input:
    path(vcf)
    path(filter_flags_variant)
    path(filter_flags_invariant)

    output:
    path("${vcf.simpleName}.vcf.gz")

    script:
    """
    echo "--gzvcf ${vcf}" > variant.args
    echo "--remove-indels" >> variant.args
    echo "--recode-INFO-all" >> variant.args
    echo "--recode" >> variant.args
    echo "--stdout" >> variant.args
    echo "--mac 1" >> variant.args
    cat ${filter_flags_variant} >> variant.args

    echo "--gzvcf ${vcf}" > invariant.args
    echo "--remove-indels" >> invariant.args
    echo "--recode-INFO-all" >> invariant.args
    echo "--recode" >> invariant.args
    echo "--stdout" >> invariant.args
    echo "--max-maf 0" >> invariant.args
    cat ${filter_flags_invariant} >> invariant.args

    cat variant.args \
        | xargs vcftools \
        | bcftools view \
            --threads ${task.cpus} \
            --exclude 'ALT="*" || TYPE!="snp"' \
            --output-type z --output ${vcf.simpleName}_variant.vcf.gz
    
    cat invariant.args \
        | xargs vcftools \
        | bcftools view \
            --threads ${task.cpus} \
            --exclude 'TYPE!="ref"' \
            --output-type z --output ${vcf.simpleName}_invariant.vcf.gz
    
    bcftools index ${vcf.simpleName}_variant.vcf.gz
    bcftools index ${vcf.simpleName}_invariant.vcf.gz

    bcftools concat \
        --allow-overlaps \
        --output-type z \
        --output ${vcf.simpleName}_tmp.vcf.gz \
        ${vcf.simpleName}_variant.vcf.gz ${vcf.simpleName}_invariant.vcf.gz

    mv ${vcf.simpleName}_tmp.vcf.gz ${vcf.simpleName}.vcf.gz
    """
}
