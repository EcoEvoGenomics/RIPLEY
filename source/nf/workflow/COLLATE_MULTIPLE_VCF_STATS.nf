include { METADATA_PREPEND_KEY_COLUMN_WITH_HEADER as PREPEND_POP_COLUMN } from "../process/metadata.nf"
include { PLOT_VCFTOOLS_VCF_STATS_POPWISE } from "../process/plot.nf"

workflow COLLATE_MULTIPLE_VCF_STATS {

    take:
    vcf_stats
    population_metadata

    main:
    with_key = vcf_stats
        .flatten()
        .map { it ->
            def key = it.simpleName.tokenize("_")[1]
            tuple("POP", key, it)
        } \
        | PREPEND_POP_COLUMN
    
    across_keys = with_key
        .flatten()
        .map { it ->
            def ext = it.extension
            def key = it.simpleName.tokenize("_")[0]
            tuple(key, ext, it)
        }
        .collectFile( { it -> ["${it[0]}_popwise.${it[1]}", it[2]] },
            skip: 1,
            keepHeader: true
        )
        .collect()

    plot = PLOT_VCFTOOLS_VCF_STATS_POPWISE(across_keys, population_metadata)

    emit:
    data = across_keys
    plot = plot

}
