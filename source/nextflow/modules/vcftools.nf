process VCFTOOLS_EXCLUDE_BED {

    label "VCFTOOLS"

    input:
    path(vcf)
    path(bed)

    output:
    path("${vcf.simpleName}_bed_excluded.vcf.gz")

    script:
    """
    vcftools --gzvcf ${vcf} --exclude-bed ${bed} --recode --stdout | gzip -c > ${vcf.simpleName}_bed_excluded.vcf.gz
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
        path("${vcf.simpleName}.ldepth.mean"), \
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
    """
}

process VCFTOOLS_CALCULATE_PAIRWISE_FST {

    label "VCFTOOLS"

    input:
    tuple path(vcf), path(pop_a_list), path(pop_b_list)

    output:
    path("${pop_a_list.simpleName}_${pop_b_list.simpleName}.weir.fst"), emit: full
    path("${pop_a_list.simpleName}_${pop_b_list.simpleName}.out"), emit: means

    script:
    """
    vcftools --gzvcf "${vcf}" \
    --weir-fst-pop "${pop_a_list}" \
    --weir-fst-pop "${pop_b_list}" \
    --out "./${pop_a_list.simpleName}_${pop_b_list.simpleName}" \
    &> "./${pop_a_list.simpleName}_${pop_b_list.simpleName}.out"
    """
}

process VCFTOOLS_CALCULATE_RELATEDNESS {

    label "VCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.relatedness2"), emit: relatedness

    script:
    """
    vcftools --gzvcf "${vcf}" \
    --relatedness2 \
    --out "./${vcf.simpleName}"
    """
}
