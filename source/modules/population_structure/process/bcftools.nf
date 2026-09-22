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
