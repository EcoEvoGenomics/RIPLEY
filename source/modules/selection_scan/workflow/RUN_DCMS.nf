include { CALCULATE_DCMS } from "../process/dcms.nf"
include { PLOT_DCMS } from "../process/plot.nf"

workflow RUN_DCMS {

    take:
    dcms_formatted
    population_metadata
    chrom_labels

    main:
    // Guard sits on the dataflow path so CALCULATE_DCMS cannot start before it lands
    collated = dcms_formatted
        .collect()
        .map { tables ->
            def stats = tables.collect { table -> table.name.tokenize(".")[-3] }.unique()
            if (stats.size() < 2) {
                error("DCMS requires at least two statistics to composite. Only ${stats.join(', ')} was formatted for it.")
            }
            return tables
        }

    dcms = CALCULATE_DCMS(file("${moduleDir}/../library/calculate_dcms.R", checkIfExists: true), collated)
    tables = dcms.collect()
    plot = PLOT_DCMS(
        file("${moduleDir}/../library/plot_dcms_main.R", checkIfExists: true),
        tables,
        population_metadata,
        chrom_labels
    )

    emit:
    data = dcms
    plot = plot

}
