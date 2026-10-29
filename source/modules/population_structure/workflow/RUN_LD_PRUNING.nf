include { PLINK_LD_PRUNE; PLINK_EXTRACT_SITES; PLINK_FILTER_MISSING_SAMPLES } from "../process/plink.nf"

workflow RUN_LD_PRUNING {

    take:
    plinkfiles
    window_kb
    step_snps
    threshold
    mind

    main:
    channel.of(mind).subscribe { m ->
        (m as Double) > 0 && (m as Double) <= 1 ?: error("prune_mind_threshold must be greater than 0 and at most 1.")
    }

    pruned = PLINK_LD_PRUNE(plinkfiles, window_kb, step_snps, threshold)
    keep = pruned.pruned_in
    extracted = PLINK_EXTRACT_SITES(plinkfiles, keep)

    // Missingness is measured on the pruned sites only
    filtered = PLINK_FILTER_MISSING_SAMPLES(extracted, mind)

    emit:
    plinkfiles = filtered.plinkfiles
    removed = filtered.removed
    logs = filtered.log

}
