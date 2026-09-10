include { BCFTOOLS_SELECT_CHROMS } from "../process/bcftools.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../process/vcftools.nf"
include { PLINK_INIT_PLINKFILES; PLINK_TO_VCF } from "../process/plink.nf"

workflow PARSE_VCF {

    // To-do: Add test for strictly alphanumeric input VCF names (e.g. genotypes.vcf.gz, chr1.vcf.gz, chr2.vcf.gz, etc.)

    take:
    vcf_path
    exclude_coords
    total_chroms
    chrom_names
    permit_solo
    permit_dir

    main:
    def input_is_solo = (file(vcf_path).isFile() && file(vcf_path).name.contains(".vcf"))
    def input_is_dir = file(vcf_path).isDirectory()
    if (input_is_solo && input_is_dir) { exit(1, "The input path may be interpreted both as file and directory.") }
    if (!(input_is_solo || input_is_dir)) { exit(1, "The input path does not exist or is not a directory or VCF.") }
    if (input_is_solo && !permit_solo) { exit(1, "This pipeline cannot process a single VCF file, only a directory.") }
    if (input_is_dir && !permit_dir) { exit(1, "This pipeline cannot process a directory, only a single VCF file.") }
    keep_chroms = chrom_names.collect()

    // Directory input assumes one vcf corresponds to exactly one chromosome
    if (input_is_dir) {
        vcf_annotated = Channel.fromPath("${vcf_path}/**.vcf.gz")
            .combine(keep_chroms.toList())
            .filter { i ->
                def vcf = i[0]
                def chroms = i[1]
                chroms.any { chrom -> vcf.simpleName.contains(chrom) }
            }
            .map { i -> i[0] }
            .ifEmpty { exit(1, "Path ${vcf_path} contains no vcf.gz files.") }
    }

    if (input_is_solo) {
        chrom_flag = keep_chroms.map { i -> i.join(",") }
        vcf_annotated = BCFTOOLS_SELECT_CHROMS(vcf_path, chrom_flag)
    }

    def exclude_path = exclude_coords ?: null
    if (exclude_path != null) {
        vcf_filtered = VCFTOOLS_EXCLUDE_BED(vcf_annotated, file(exclude_path, checkIfExists: true))
    } else {
        vcf_filtered = vcf_annotated
    }

    plinkfiles = PLINK_INIT_PLINKFILES(vcf_filtered, total_chroms)
    vcf_condensed = PLINK_TO_VCF(plinkfiles)

    emit:
    as_plinkfiles = plinkfiles
    vcf_condensed = vcf_condensed   // PLINK condenses VCFs by removing annotations
    vcf_annotated = vcf_filtered    // ... but the annotations are sometimes useful

}
