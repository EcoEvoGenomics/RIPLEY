include { PLINK_LD_PRUNE; PLINK_EXTRACT_SITES } from "../modules/plink.nf"

workflow RUN_LD_PRUNING {

    take:
    bedfiles
    window_kb
    step_snps
    threshold

    main:
    pruned = PLINK_LD_PRUNE(bedfiles, window_kb, step_snps, threshold)
    keep = pruned.prune_in
    drop = pruned.prune_out
    bedfiles_out = PLINK_EXTRACT_SITES(bedfiles, keep)

    emit:
    keep = keep
    drop = drop
    bedfiles = bedfiles_out

}
