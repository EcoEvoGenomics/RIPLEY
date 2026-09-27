include { VCFTOOLS_FILTER_VARIANTS } from "../process/vcftools.nf"
include { PARSE_FILTER_FLAGS } from "./PARSE_FILTER_FLAGS.nf"

workflow RUN_VCF_FILTERING {

    take:
    vcfs
    filter_flags_path

    main:
    filter_flags = PARSE_FILTER_FLAGS(filter_flags_path).flags
    filtered = VCFTOOLS_FILTER_VARIANTS(vcfs, filter_flags)

    emit:
    vcf_filtered = filtered
    flags = filter_flags

}
