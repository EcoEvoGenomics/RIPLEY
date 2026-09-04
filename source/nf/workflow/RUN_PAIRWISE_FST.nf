include { PAIR_CHANNEL_TO_SELF } from "./PAIR_CHANNEL_TO_SELF.nf"
include { WRITE_POPULATION_CENSUS } from "../process/system.nf"
include { VCFTOOLS_CALCULATE_PAIRWISE_FST } from "../process/vcftools.nf"
include { PLOT_VCFTOOLS_PAIRWISE_MEAN_FST } from "../process/plotting.nf"

workflow RUN_PAIRWISE_FST {

    take:
    vcf
    pop_list
    metadata

    main:
    pop_censuses = WRITE_POPULATION_CENSUS(pop_list, metadata)
    pairwise_pop_censuses = PAIR_CHANNEL_TO_SELF(pop_censuses)
    results = VCFTOOLS_CALCULATE_PAIRWISE_FST(vcf.combine(pairwise_pop_censuses))
    mean = results.mean
        .map { result -> 
            def pop_a = result[0] as String
            def pop_b = result[1] as String
            def fst   = result[2] as Double
            if (fst <= 0) { fst = 0 }
            "${pop_a}\t${pop_b}\t${fst}\n"
        }
        .collectFile( name: "weighted.fst", sort: { pop_pair -> pop_pair[0] } )
    plot = PLOT_VCFTOOLS_PAIRWISE_MEAN_FST(mean)

    emit:
    logfile = results.logfile
    data = results.full
    mean = mean
    plot = plot

}
