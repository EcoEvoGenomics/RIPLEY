include { SPLIT_VCF } from "./SPLIT_VCF.nf"
include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { DROP_MISMATCHED_FILEKEY_PAIRS } from "./DROP_MISMATCHED_FILEKEY_PAIRS.nf"
include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../process/genomicsgeneral.nf"
include { BCFTOOLS_MERGE_VCFS; BCFTOOLS_LIST_SAMPLES } from "../process/bcftools.nf"

workflow RUN_POPGEN_WINDOWS {

    take:
    vcfs
    metadata
    window_size
    step_size
    min_sites

    main:
    repo = GET_GENOMICS_GENERAL()

    pairwise_vcf = PAIR_CHANNEL_TO_SELF(vcfs) | DROP_MISMATCHED_FILEKEY_PAIRS
    merged_vcf = BCFTOOLS_MERGE_VCFS(pairwise_vcf)
    merged_geno = GENOMICS_GENERAL_VCF_TO_GENO(repo, merged_vcf)
    target_samples = BCFTOOLS_LIST_SAMPLES(merged_vcf)

    inputs = merged_geno
        .combine(target_samples)
        .map { files ->
            def key_geno = files[0].simpleName
            def key_samples = files[1].simpleName
            key_geno == key_samples
                ? tuple(files[0], files[1])
                : null
        }
        .combine(metadata)

    GENOMICS_GENERAL_POPGEN_WINDOWS(repo, inputs, window_size, step_size, min_sites)

    emit:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
    
}
