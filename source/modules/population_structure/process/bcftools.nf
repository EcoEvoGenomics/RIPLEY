process BCFTOOLS_VCF_TO_GENOTABLE {

    label "BCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.gt")

    script:
    """
    bcftools query -l "${vcf}" | paste -sd'\\t' | sed 's/^/ID\\t/' > "${vcf.simpleName}.gt"
    bcftools query -f '%ID[\\t%TGT]\\n' "${vcf}" >> "${vcf.simpleName}.gt"
    """
}

process BCFTOOLS_LIST_SITE_IDS {

    // PLINK may rename contigs on VCF export, so CHROM:POS is not a safe proxy for the variant ID

    label "BCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.sitekeys")

    script:
    """
    bcftools query -f '%CHROM\\t%POS\\t%ID\\n' "${vcf}" > "${vcf.simpleName}.sitekeys"
    """
}
