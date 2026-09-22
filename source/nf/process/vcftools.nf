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

process VCFTOOLS_SNP_DENSITY {

    label "VCFTOOLS"

    input:
    path(vcf)
    val(binsize)

    output:
    path("${vcf.simpleName}.snpden")

    script:
    """
    vcftools --gzvcf ${vcf} --SNPdensity ${binsize} --out ${vcf.simpleName}
    sed -i -E 's/^([0-9]+)\\t/chr\\1\\t/' ${vcf.simpleName}.snpden
    """
}

process VCFTOOLS_VCF_STATS {

    label "VCFTOOLS"

    input:
    path(vcf)

    output:
    tuple \
        path("${vcf.simpleName}.frq"), \
        path("${vcf.simpleName}.idepth"), \
        path("${vcf.simpleName}.imiss"), \
        path("${vcf.simpleName}.ldepth"), \
        path("${vcf.simpleName}.lqual"), \
        path("${vcf.simpleName}.lmiss"), \
        path("${vcf.simpleName}.het"), \
        path("${vcf.simpleName}.hwe")

    script:
    """
    vcftools --gzvcf ${vcf} --freq2 --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --depth --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --missing-indv --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --site-mean-depth --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --site-quality --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --missing-site --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --het --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --hardy --out ${vcf.simpleName}
    mv ${vcf.simpleName}.ldepth.mean ${vcf.simpleName}.ldepth
    """
}

