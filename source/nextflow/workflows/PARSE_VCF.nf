include { BCFTOOLS_SELECT_CHROMS } from "../modules/bcftools.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../modules/vcftools.nf"
include { PLINK_INIT_BEDFILES; PLINK_TO_VCF } from "../modules/plink.nf"

workflow PARSE_VCF {

    take:
    vcf_path
    exclude_coords
    n_chroms
    chroms
    permit_dir

    main:

    if (file(vcf_path).isDirectory()) {

        if (!permit_dir) {
            exit(1, "This pipeline can only process a single VCF file.")
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

    def exclude_path = exclude_coords ?: null

    if (exclude_path == null) {
        vcf_filtered = vcf_annotated
    } else {
        vcf_filtered = VCFTOOLS_EXCLUDE_BED(vcf_annotated, file(exclude_path, checkIfExists: true))
    } //

    bedfiles = PLINK_INIT_BEDFILES(vcf_filtered, plink_n_chroms)
    vcf_condensed = PLINK_TO_VCF(bedfiles)

    emit:
    bedfiles = bedfiles
    vcf = vcf_condensed
    vcf_annotated = vcf_annotated
    vcf = vcf_condensed             // PLINK condenses VCFs by removing annotations
    vcf_annotated = vcf_filtered    // ... but the annotations are sometimes useful

}
