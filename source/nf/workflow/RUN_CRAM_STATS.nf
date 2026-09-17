include { SAMTOOLS_STATS } from "../process/samtools.nf"
include { METADATA_PREPEND_KEY_COLUMN as PREPEND_SAMPLE_COLUMN } from "../process/metadata.nf"
include { PLOT_SAMTOOLS_CRAM_STATS } from "../process/plot.nf"

workflow RUN_CRAM_STATS {

    take:
    cram_indexed

    main:
    stats = SAMTOOLS_STATS(cram_indexed)

    stats_with_key = stats
        .flatten()
        .map { it ->
            def key = it.simpleName
            tuple(key, it)
        } | PREPEND_SAMPLE_COLUMN

    stats_combined = stats_with_key
        .map { it ->
            def ext = it.extension
            tuple(ext, it)
        }
        .collectFile( { it -> ["stats.${it[0]}", it[1]]} )

    stats_plot = PLOT_SAMTOOLS_CRAM_STATS(stats_combined)

    emit:
    data = stats_combined
    plot = stats_plot

}
