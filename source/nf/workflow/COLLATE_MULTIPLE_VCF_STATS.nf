include { METADATA_PREPEND_KEY_COLUMN_WITH_HEADER } from "../process/metadata.nf"
include { PLOT_COLLATED_VCF_STATS } from "../process/plot.nf"

workflow COLLATE_MULTIPLE_VCF_STATS {

    take:
    vcf_stats

    main:
    with_key = vcf_stats
        .flatten()
        .map { it ->
            def key = it.simpleName.tokenize("_")[1]
            tuple("KEY", key, it)
        } \
        | METADATA_PREPEND_KEY_COLUMN_WITH_HEADER
    
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

    plot = PLOT_COLLATED_VCF_STATS(across_keys)

    emit:
    data = across_keys
    plot = plot

}
