include { PLINK_LD_PRUNE; PLINK_EXTRACT_SITES } from "../modules/plink.nf"

workflow RUN_LD_PRUNING {

    take:
    plinkfiles
    window_kb
    step_snps
    threshold

    main:
    pruned = PLINK_LD_PRUNE(plinkfiles, window_kb, step_snps, threshold)
    keep = pruned.pruned_in
    drop = pruned.pruned_out
    plinkfiles_out = PLINK_EXTRACT_SITES(plinkfiles, keep)

    emit:
    keep = keep
    drop = drop
    plinkfiles = plinkfiles_out

}
