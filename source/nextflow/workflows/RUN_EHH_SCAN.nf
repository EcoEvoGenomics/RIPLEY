include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { WRITE_POPULATION_CENSUS } from "../modules/system.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../modules/vcftools.nf"
include { BCFTOOLS_PICK_SAMPLES } from "../modules/bcftools.nf"
include { REHH_LOAD_VCF; REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY } from "../modules/rehh.nf"
include { REHH_CALCULATE_IHS; REHH_CALCULATE_XPEHH } from "../modules/rehh.nf"
include { PLOT_REHH_XPEHH } from "../modules/plotting.nf"

workflow RUN_EHH_SCAN {

    take:
    vcfs
    metadata
    populations
    window_size
    step_size
    min_sites

    main:
    population_vcfs = WRITE_POPULATION_CENSUS(populations, metadata) | combine(vcfs) | BCFTOOLS_PICK_SAMPLES
    population_ehh = REHH_LOAD_VCF(population_vcfs) | REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY | flatten

    pairwise_ehh = PAIR_CHANNEL_TO_SELF(population_ehh)
        .map { pair ->
            def key_a = file(pair[0]).SimpleName.toString().tokenize("_")[0]
            def key_b = file(pair[1]).SimpleName.toString().tokenize("_")[0]
            def population_a = file(pair[0]).SimpleName.toString().tokenize("_")[1]
            def population_b = file(pair[1]).SimpleName.toString().tokenize("_")[1]
            key_a == key_b
                ? tuple(key_a, population_a, population_b, file(pair[0]), file(pair[1]))
                : null // Do not compare across keys (keys correspond to original vcfs)
        }
    
    REHH_CALCULATE_IHS(population_ehh, window_size, step_size, min_sites)
    REHH_CALCULATE_XPEHH(pairwise_ehh, window_size, step_size, min_sites)

    emit:
    ihs = REHH_CALCULATE_IHS.out
    xpehh = REHH_CALCULATE_XPEHH.out
    
}
