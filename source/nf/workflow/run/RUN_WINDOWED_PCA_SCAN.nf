include { PAIR_CHANNEL_TO_SELF } from "../utils/PAIR_CHANNEL_TO_SELF.nf"
include { DROP_MISMATCHED_FILEKEY_PAIRS } from "../utils/DROP_MISMATCHED_FILEKEY_PAIRS.nf"
include { BCFTOOLS_MERGE_VCFS } from "../../process/bcftools.nf"
include { GET_WINPCA; WINPCA_PCA_CHROMWISE } from "../../process/winpca.nf"

workflow RUN_WINDOWED_PCA_SCAN {

    take:
    vcfs
    chrom_indices
    window_size
    step_size

    main:
    repo = GET_WINPCA()

    pairwise_vcf = PAIR_CHANNEL_TO_SELF(vcfs) | DROP_MISMATCHED_FILEKEY_PAIRS
    merged_vcf = pairwise_vcf
        .map { key, name_a, name_b, vcf_a, vcf_b -> tuple("${key}_${name_a}_${name_b}", [vcf_a, vcf_b]) }
        | BCFTOOLS_MERGE_VCFS

    wpca_inputs = merged_vcf.combine(chrom_indices)
        .filter { it -> it[0].simpleName.tokenize("_")[0] == it[1] }
        .map { it ->
            def vcf = it[0]
            def chrom = it[1]
            def chrom_length = it[2]
            tuple (vcf, chrom, chrom_length, window_size, step_size)
        }

    wpca = WINPCA_PCA_CHROMWISE(repo, wpca_inputs)

    emit:
    data = wpca
    
}
