include { PAIR_CHANNEL_TO_SELF } from "../utils/PAIR_CHANNEL_TO_SELF.nf"
include { DROP_MISMATCHED_FILEKEY_PAIRS } from "../utils/DROP_MISMATCHED_FILEKEY_PAIRS.nf"
include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../../process/genomics.nf"
include { BCFTOOLS_MERGE_VCFS; BCFTOOLS_LIST_SAMPLES } from "../../process/bcftools.nf"

workflow RUN_POPGEN_WINDOWS_SCAN {

    take:
    vcfs
    sample_metadata
    window_size
    step_size
    min_sites

    main:
    repo = GET_GENOMICS_GENERAL()

    pairwise_vcf = PAIR_CHANNEL_TO_SELF(vcfs) | DROP_MISMATCHED_FILEKEY_PAIRS
    merged_vcf = pairwise_vcf
        .map { key, name_a, name_b, vcf_a, vcf_b -> tuple("${key}_${name_a}_${name_b}", [vcf_a, vcf_b]) }
        | BCFTOOLS_MERGE_VCFS
    merged_geno = GENOMICS_GENERAL_VCF_TO_GENO(repo, merged_vcf)
    target_samples = BCFTOOLS_LIST_SAMPLES(merged_vcf)

    inputs = merged_geno
        .combine(target_samples)
        .filter { files -> files[0].simpleName == files[1].simpleName }
        .map { files -> tuple(files[0], files[1]) }
        .combine(sample_metadata)

    GENOMICS_GENERAL_POPGEN_WINDOWS(repo, inputs, window_size, step_size, min_sites)

    emit:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
    
}
