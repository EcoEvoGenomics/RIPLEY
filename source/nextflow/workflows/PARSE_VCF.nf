include { BCFTOOLS_SELECT_CHROMS } from "../modules/bcftools.nf"
include { PLINK_INIT_BEDFILES; PLINK_TO_VCF } from "../modules/plink.nf"

workflow PARSE_VCF {

    take:
    vcf_path
    n_chroms
    chroms
    permit_dir

    main:

    if (file(vcf_path).isDirectory()) {

        if (!permit_dir) {
            exit(1, "ERROR: This pipeline can only process a single VCF file.")
        }

        vcf_annotated = Channel.fromPath("${vcf_path}/**.vcf.gz")
            .combine(chroms)
            .filter { i -> i[1].tokenize(",").any { j -> i[0].simpleName.contains(j) } }  // For any VCF of chrom "X" assume name uniquely contains "X"
            .map { i -> i[0] }
        bedfiles = PLINK_INIT_BEDFILES(vcf_annotated, 1)
        vcf_condensed = PLINK_TO_VCF(bedfiles)

    } else if (file(vcf_path).isFile()) {

        vcf_annotated = BCFTOOLS_SELECT_CHROMS(vcf_path, chroms)
        bedfiles = PLINK_INIT_BEDFILES(vcf_annotated, n_chroms)
        vcf_condensed = PLINK_TO_VCF(bedfiles)

    }

    emit:
    bedfiles = bedfiles
    vcf = vcf_condensed
    vcf_annotated = vcf_annotated

}
