include { BCFTOOLS_EXCLUDE_CHROMS } from "../modules/bcftools.nf"
include { PLINK_INIT_BEDFILES; PLINK_TO_VCF } from "../modules/plink.nf"

workflow PARSE_VCF_FILE {

    take:
    vcf_path
    n_chroms
    exclude_chroms

    main:
    vcf_annotated = BCFTOOLS_EXCLUDE_CHROMS(vcf_path, exclude_chroms)
    bedfiles = PLINK_INIT_BEDFILES(vcf_annotated, n_chroms)
    vcf_condensed = PLINK_TO_VCF(bedfiles)

    emit:
    vcf = vcf_condensed
    vcf_annot = vcf_annotated
    bedfiles = bedfiles
}
