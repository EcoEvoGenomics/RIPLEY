include { VCFTOOLS_FILTER_VARIANTS_AND_INVARIANTS } from "../process/vcftools.nf"
include { PARSE_FILTER_FLAGS as PARSE_VARIANT_FLAGS; PARSE_FILTER_FLAGS as PARSE_INVARIANT_FLAGS } from "./PARSE_FILTER_FLAGS.nf"

workflow RUN_VCF_FILTERING_VARIANT_INVARIANT {

    take:
    vcf
    flags_variant_path
    flags_invariant_path

    main:
    flags_variant = PARSE_VARIANT_FLAGS(flags_variant_path).flags
    flags_invariant = PARSE_INVARIANT_FLAGS(flags_invariant_path).flags
    filtered = VCFTOOLS_FILTER_VARIANTS_AND_INVARIANTS(vcf, flags_variant, flags_invariant)

    emit:
    vcf_filtered = filtered
    flags = flags_variant.mix(flags_invariant)

}
