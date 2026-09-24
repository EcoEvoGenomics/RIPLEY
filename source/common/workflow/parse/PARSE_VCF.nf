include { BCFTOOLS_FILTER_CHROMS; BCFTOOLS_INDEX; BCFTOOLS_COUNT_RECORDS } from "../../process/bcftools.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../../process/vcftools.nf"
include { alphanumericIssue; keyFor } from "../../library/filekeys.nf"

workflow PARSE_VCF {

    take:
    vcf_path
    exclude_coords
    chrom_names
    permit_solo
    permit_dir

    main:
    def permitted_extensions = ["vcf.gz", "vcf"]

    def input_file = file(vcf_path)
    if (!(input_file instanceof Path)) { error("The input path must name a single VCF or directory, not a glob pattern: ${vcf_path}") }

    def vcf_name = input_file.name
    def input_is_solo = (input_file.isFile() && keyFor(vcf_name, permitted_extensions) != null)
    def input_is_dir = input_file.isDirectory()
    if (input_is_solo && input_is_dir) { error("The input path may be interpreted both as file and directory: ${vcf_path}") }
    if (!(input_is_solo || input_is_dir)) { error("The input path does not exist or is not a directory or VCF: ${vcf_path}") }
    if (input_is_solo && !permit_solo) { error("This pipeline cannot process a single VCF file, only a directory: ${vcf_path}") }
    if (input_is_dir && !permit_dir) { error("This pipeline cannot process a directory, only a single VCF file: ${vcf_path}") }
    keep_chroms = chrom_names.collect()

    // Directory input expects one vcf per chromosome (chr1.vcf.gz, chr2.vcf.gz, ... chrN.vcf.gz), verified against contents below
    if (input_is_dir) {
        vcf_found = Channel.fromPath("${vcf_path}/**.vcf.gz")
            .ifEmpty { error("Path ${vcf_path} contains no vcf.gz files.") }

        // The whole name minus its extension is the chromosome key: downstream tokenising requires strict alphanumeric names
        vcf_found.subscribe { vcf ->
            def issue = alphanumericIssue(keyFor(vcf.name, permitted_extensions), "chromosome key", "Input VCF ${vcf.name}")
            if (issue) { error(issue) }
        }

        vcf_retained = vcf_found
            .combine(keep_chroms.toList())
            .filter { i ->
                def vcf = i[0]
                def chroms = i[1]
                chroms.contains(keyFor(vcf.name, permitted_extensions))
            }
            .map { i -> i[0] }
    }

    if (input_is_solo) {
        def issue = alphanumericIssue(keyFor(vcf_name, permitted_extensions), "file key", "Input VCF ${vcf_name}")
        if (issue) { error(issue) }
        chrom_flag = keep_chroms.map { i -> i.join(",") }
        vcf_retained = BCFTOOLS_FILTER_CHROMS(vcf_path, chrom_flag)
    }

    def exclude_path = exclude_coords ?: null
    if (exclude_path != null) {
        vcf_filtered = VCFTOOLS_EXCLUDE_BED(vcf_retained, file(exclude_path, checkIfExists: true))
    } else {
        vcf_filtered = vcf_retained
    }

    vcf_indexed = BCFTOOLS_INDEX(vcf_filtered)
    vcf_nrecords = BCFTOOLS_COUNT_RECORDS(vcf_indexed).per_chrom
        .splitCsv(sep: "\t")
        .filter { row -> row[2].toInteger() > 0 }
        .map { row -> row[0] }

    // Empty or misnamed VCFs fail here
    vcf_nrecords
        .collect()
        .ifEmpty([])
        .map { found -> [found] }
        .combine(keep_chroms.toList())
        .subscribe { i ->
            def missing = i[1] - i[0]
            if (missing) {
                error("Input ${vcf_path} has no variants on retained chromosome(s): ${missing.join(', ')}.")
            }
        }

    emit:
    vcf = vcf_filtered
    vcf_indexed = vcf_indexed

}
