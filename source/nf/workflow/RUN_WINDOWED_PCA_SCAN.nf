include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { DROP_MISMATCHED_FILEKEY_PAIRS } from "./DROP_MISMATCHED_FILEKEY_PAIRS.nf"
include { BCFTOOLS_MERGE_VCFS } from "../process/bcftools.nf"
include { GET_WINPCA; WINPCA_PCA_CHROMWISE } from "../process/winpca.nf"

workflow RUN_WINDOWED_PCA_SCAN {

    take:
    vcfs
    chrom_indices
    window_size
    step_size

    main:
    repo = GET_WINPCA()

    pairwise_vcf = PAIR_CHANNEL_TO_SELF(vcfs) | DROP_MISMATCHED_FILEKEY_PAIRS
    merged_vcf = BCFTOOLS_MERGE_VCFS(pairwise_vcf)

    wpca = WINPCA_PCA_CHROMWISE(
        repo,
        merged_vcf.combine(chrom_indices)
            .map { it -> 
                def vcf = it[0]
                def vcf_key = vcf.simpleName.tokenize("_")[0]
                def chrom = it[1]
                def chrom_length = it[2]
                vcf_key == chrom
                    ? tuple (vcf, chrom, chrom_length, window_size, step_size)
                    : null
            }
    )

    emit:
    data = wpca
    
}
