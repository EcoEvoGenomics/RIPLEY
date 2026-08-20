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
    
    def input_is_vcf = (file(vcf_path).isFile() && file(vcf_path).name.contains(".vcf"))
    def input_is_dir = file(vcf_path).isDirectory()

    if (input_is_vcf && input_is_dir) {
        exit(1, "The input path may be interpreted both as file and directory.")
    }

    if (!(input_is_vcf || input_is_dir)) {
        exit(1, "The input path does not exist or is not a directory or VCF.")
    }

    if (input_is_dir && !permit_dir) {
        exit(1, "This pipeline cannot process a directory, only a single VCF file.")
    }

    if (input_is_dir) {

        // For any VCF of chrom "X" assume name uniquely contains "X"
        vcf_annotated = Channel.fromPath("${vcf_path}/**.vcf.gz", checkIfExist: true)
            .combine(chroms)
            .filter { i -> i[1].tokenize(",").any { j -> i[0].simpleName.contains(j) } }
            .map { i -> i[0] }
        
        // Assume one chromosome in each VCF
        plink_n_chroms = Channel.value(1)

    }

    if (input_is_vcf) {
        vcf_annotated = BCFTOOLS_SELECT_CHROMS(vcf_path, chroms)
        plink_n_chroms = n_chroms
    }

    def exclude_path = exclude_coords ?: null

    if (exclude_path == null) {
        vcf_filtered = vcf_annotated
    } else {
        vcf_filtered = VCFTOOLS_EXCLUDE_BED(vcf_annotated, file(exclude_path, checkIfExists: true))
    }

    bedfiles = PLINK_INIT_BEDFILES(vcf_filtered, plink_n_chroms)
    vcf_condensed = PLINK_TO_VCF(bedfiles)

    emit:
    bedfiles = bedfiles
    vcf = vcf_condensed             // PLINK condenses VCFs by removing annotations
    vcf_annotated = vcf_filtered    // ... but the annotations are sometimes useful

}
