include { PAIR_CHANNEL_TO_SELF } from "../../../common/workflow/utils/PAIR_CHANNEL_TO_SELF.nf"
include { VCFTOOLS_CALCULATE_PAIRWISE_FST } from "../../process/vcftools.nf"
include { PLOT_VCFTOOLS_PAIRWISE_MEAN_FST } from "../../process/plot.nf"

workflow RUN_PAIRWISE_FST {

    take:
    vcf
    pop_censuses

    main:
    pairwise_pop_censuses = PAIR_CHANNEL_TO_SELF(pop_censuses)
    results = VCFTOOLS_CALCULATE_PAIRWISE_FST(vcf.combine(pairwise_pop_censuses))
    
    mean = results.mean
        // Result is NA when a pair yields no usable sites, which is not castable to Double
        .filter { result ->
            def valid_result = (result[2] as String) ==~ /^-?\d+(\.\d+)?([eE][-+]?\d+)?$/
            valid_result
        }
        .map { result -> 
            def pop_a = result[0] as String
            def pop_b = result[1] as String
            def fst   = result[2] as Double
            if (fst <= 0) { fst = 0 }
            "${pop_a}\t${pop_b}\t${fst}\n"
        }
        .collectFile( name: "weighted.fst", sort: { line -> line.tokenize("\t")[0] } )
    
    plot = PLOT_VCFTOOLS_PAIRWISE_MEAN_FST(
        file("${moduleDir}/../../../R/plot_fst.R", checkIfExists: true),
        mean
    )

    emit:
    logs = results.logs
    data = results.full
    mean = mean
    plot = plot

}
