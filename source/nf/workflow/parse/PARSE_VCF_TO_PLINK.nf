include { PLINK_INIT_PLINKFILES; PLINK_TO_VCF } from "../../process/plink.nf"

workflow PARSE_VCF_TO_PLINK {

    take:
    vcf
    total_chroms

    main:
    plinkfiles = PLINK_INIT_PLINKFILES(vcf, total_chroms)
    vcf_condensed = PLINK_TO_VCF(plinkfiles)

    emit:
    as_plinkfiles = plinkfiles
    vcf_condensed = vcf_condensed   // PLINK condenses VCFs by dropping annotations (and indels and multiallelic sites)

}
