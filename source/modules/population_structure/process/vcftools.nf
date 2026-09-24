process VCFTOOLS_CALCULATE_PAIRWISE_FST {

    label "VCFTOOLS"

    input:
    tuple path(vcf), path(pop_a_list), path(pop_b_list)

    output:
    path("${pop_a_list.simpleName}_${pop_b_list.simpleName}.weir.fst"), emit: full
    path("${pop_a_list.simpleName}_${pop_b_list.simpleName}.out"), emit: logs
    tuple \
        val("${pop_a_list.simpleName}"), \
        val("${pop_b_list.simpleName}"),
        env("weighted_mean_fst"), emit: mean

    script:
    """
    vcftools --gzvcf "${vcf}" \
    --weir-fst-pop "${pop_a_list}" \
    --weir-fst-pop "${pop_b_list}" \
    --out "./${pop_a_list.simpleName}_${pop_b_list.simpleName}" \
    &> "./${pop_a_list.simpleName}_${pop_b_list.simpleName}.out"

    weighted_mean_fst=\$(grep 'Weir and Cockerham weighted Fst estimate:' \
        "./${pop_a_list.simpleName}_${pop_b_list.simpleName}.out" | awk '{print \$NF}')
    weighted_mean_fst=\${weighted_mean_fst:-NA}
    export weighted_mean_fst
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

process VCFTOOLS_CALCULATE_ALLELE_FREQUENCIES {

    // --freq rather than --freq2 because allele identity is needed downstream.

    label "VCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.frq"), emit: frequencies

    script:
    """
    vcftools --gzvcf "${vcf}" \
    --freq \
    --out "./${vcf.simpleName}"
    """
}
