include { BCFTOOLS_INDEX; BCFTOOLS_INDEX as BCFTOOLS_INDEX_INPUT; BCFTOOLS_PICK_CHROM; BCFTOOLS_COUNT_RECORDS; BCFTOOLS_EXCLUDE_BED } from "../../process/bcftools.nf"
include { alphanumericIssue; keyFor } from "../../library/filekeys.nf"

workflow PARSE_VCF {

    // Emits one VCF per chromosome whether the input names a single VCF or a directory.

    take:
    vcf_path
    exclude_coords
    chrom_names

    main:
    def permitted_extensions = ["vcf.gz", "vcf"]

    def input_file = file(vcf_path)
    if (!(input_file instanceof Path)) { error("The input path must name a single VCF or directory, not a glob pattern: ${vcf_path}") }

    def vcf_name = input_file.name
    def input_is_solo = (input_file.isFile() && keyFor(vcf_name, permitted_extensions) != null)
    def input_is_dir = input_file.isDirectory()
    if (!(input_is_solo || input_is_dir)) { error("The input path does not exist or is not a directory or VCF: ${vcf_path}") }

    // Directory input expects one vcf per chromosome (chr1.vcf.gz, ... chrN.vcf.gz)
    if (input_is_dir) {
        vcf_found = Channel.fromPath("${vcf_path}/**.vcf.gz").ifEmpty { error("Path ${vcf_path} contains no vcf.gz files.") }
    }

    // Solo input is split to match the expected directory shape
    if (input_is_solo) {
        vcf_found = BCFTOOLS_PICK_CHROM(BCFTOOLS_INDEX_INPUT(vcf_path), chrom_names)
    }

    // The whole name minus its extension is the chromosome key: downstream tokenising requires strict alphanumeric names
    vcf_found.subscribe { vcf ->
        def issue = alphanumericIssue(keyFor(vcf.name, permitted_extensions), "chromosome key", "Input VCF ${vcf.name}")
        if (issue) { error(issue) }
    }

    // Retain only chroms retained by PARSE_REFERENCE_GENOME
    keep_chroms = chrom_names.collect()
    vcf_retained = vcf_found
        .combine(keep_chroms.toList())
        .filter { i ->
            def vcf = i[0]
            def chroms = i[1]
            chroms.contains(keyFor(vcf.name, permitted_extensions))
        }
        .map { i -> i[0] }

    // Index before exclusion to check bedfile against contigs from csi rather than from vcf
    def exclude_bed = exclude_coords ? Channel.value(file(exclude_coords, checkIfExists: true)) : null
    vcf_indexed = BCFTOOLS_INDEX(vcf_retained)
    if (exclude_bed == null) {
        vcf_indexed_out = vcf_indexed
    } else {
        vcf_indexed_out = BCFTOOLS_EXCLUDE_BED(vcf_indexed.combine(exclude_bed))
    }

    vcf_nrecords = BCFTOOLS_COUNT_RECORDS(vcf_indexed_out).per_chrom
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
    vcf = vcf_indexed_out.map { i -> i[0] }
    vcf_indexed = vcf_indexed_out

}
