include {PLINK_INIT_BEDFILES; PLINK_EXCLUDE_CHROMS; PLINK_TO_VCF} from "../modules/plink.nf"

workflow PARSE_VCF_FILE {

    take:
    vcf_path
    n_chroms
    exclude_chroms

    main:
    init_bedfiles = PLINK_INIT_BEDFILES(vcf_path, n_chroms)
    filtered_bedfiles = PLINK_EXCLUDE_CHROMS(init_bedfiles, exclude_chroms)
    vcf = PLINK_TO_VCF(filtered_bedfiles)

    emit:
    vcf = vcf
    bedfiles = filtered_bedfiles
}
