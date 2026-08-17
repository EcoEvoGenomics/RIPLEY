include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { WRITE_POPULATION_CENSUS } from "../modules/system.nf"
include { VCFTOOLS_CALCULATE_PAIRWISE_FST } from "../modules/vcftools.nf"

workflow RUN_PAIRWISE_FST {

    take:
    vcf
    pop_list
    metadata

    main:
    pop_censuses = WRITE_POPULATION_CENSUS(Channel.from(pop_list), metadata)
    pairwise_pop_censuses = PAIR_CHANNEL_TO_SELF(pop_censuses)
    results = VCFTOOLS_CALCULATE_PAIRWISE_FST(vcf.combine(pairwise_pop_censuses))

    emit:
    mean_fst = results.means
    site_fst = results.full
}
