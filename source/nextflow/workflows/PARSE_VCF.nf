include { BCFTOOLS_SELECT_CHROMS } from "../modules/bcftools.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../modules/vcftools.nf"
include { PLINK_INIT_PLINKFILES; PLINK_TO_VCF } from "../modules/plink.nf"

workflow PARSE_VCF {

    take:
    vcf_path
    exclude_coords
    chrom_names
    permit_dir

    main:
    def input_is_vcf = (file(vcf_path).isFile() && file(vcf_path).name.contains(".vcf"))
    def input_is_dir = file(vcf_path).isDirectory()
    if (input_is_vcf && input_is_dir) { exit(1, "The input path may be interpreted both as file and directory.") }
    if (!(input_is_vcf || input_is_dir)) { exit(1, "The input path does not exist or is not a directory or VCF.") }
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
        plink_n_chroms = Channel.value(1)
    }

    if (input_is_vcf) {
        chrom_flag = keep_chroms.map { i -> i.join(",")}
        vcf_annotated = BCFTOOLS_SELECT_CHROMS(vcf_path, chrom_flag)
        plink_n_chroms = chrom_names.count()
    }

    // To-do: Add test for strictly alphanumeric input VCF names (e.g. myVars.vcf.gz)

    def exclude_path = exclude_coords ?: null
    if (exclude_path != null) {
        vcf_filtered = VCFTOOLS_EXCLUDE_BED(vcf_annotated, file(exclude_path, checkIfExists: true))
    } else {
        vcf_filtered = vcf_annotated
    }

    plinkfiles = PLINK_INIT_PLINKFILES(vcf_filtered, plink_n_chroms)
    vcf_condensed = PLINK_TO_VCF(plinkfiles)

    emit:
    as_plinkfiles = plinkfiles
    vcf_condensed = vcf_condensed   // PLINK condenses VCFs by removing annotations
    vcf_annotated = vcf_filtered    // ... but the annotations are sometimes useful

}
