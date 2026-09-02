include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { DROP_MISMATCHED_FILEKEY_PAIRS } from "./DROP_MISMATCHED_FILEKEY_PAIRS.nf"
include { REHH_LOAD_VCF; REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY } from "../process/rehh.nf"
include { REHH_CALCULATE_IHS; REHH_CALCULATE_XPEHH } from "../process/rehh.nf"
include { PLOT_REHH_XPEHH } from "../process/plotting.nf"

workflow RUN_EHH_SCAN {

    take:
    vcfs
    window_size
    step_size
    min_sites

    main:
    population_ehh = REHH_LOAD_VCF(vcfs) | REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY | flatten
    pairwise_ehh = PAIR_CHANNEL_TO_SELF(population_ehh) | DROP_MISMATCHED_FILEKEY_PAIRS

    REHH_CALCULATE_IHS(population_ehh, window_size, step_size, min_sites)
    REHH_CALCULATE_XPEHH(pairwise_ehh, window_size, step_size, min_sites)

    emit:
    ihs = REHH_CALCULATE_IHS.out
    xpehh = REHH_CALCULATE_XPEHH.out
    
}
