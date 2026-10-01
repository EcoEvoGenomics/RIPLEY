include { CALCULATE_DCMS } from "../process/dcms.nf"

workflow RUN_DCMS {

    take:
    dcms_formatted

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

    emit:
    dcms = dcms

}
