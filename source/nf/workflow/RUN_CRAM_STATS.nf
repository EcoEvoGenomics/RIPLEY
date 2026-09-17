include { SAMTOOLS_STATS } from "../process/samtools.nf"
include { METADATA_PREPEND_KEY_COLUMN as PREPEND_SAMPLE_COLUMN } from "../process/metadata.nf"
include { METADATA_PREPEND_KEY_COLUMN as PREPEND_POPULATION_COLUMN } from "../process/metadata.nf"
include { PLOT_SAMTOOLS_CRAM_STATS } from "../process/plot.nf"
include { PLOT_SAMTOOLS_CRAM_STATS_POPWISE } from "../process/plot.nf"

workflow RUN_CRAM_STATS {

    take:
    cram_indexed
    population_map

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

    // Samples absent from the map are non-focal and drop out of the join
    stats_with_population = stats_with_key
        .map { it -> tuple(it.simpleName, it) }
        .join(population_map)
        .map { _sample, table, population -> tuple(population, table) } \
        | PREPEND_POPULATION_COLUMN

    stats_popwise = stats_with_population
        .map { it ->
            def ext = it.extension
            tuple(ext, it)
        }
        .collectFile( { it -> ["stats_popwise.${it[0]}", it[1]]} )

    stats_plot = PLOT_SAMTOOLS_CRAM_STATS(stats_combined)
    stats_popwise_plot = PLOT_SAMTOOLS_CRAM_STATS_POPWISE(stats_popwise)

    emit:
    data = stats_combined
    plot = stats_plot
    popwise_data = stats_popwise
    popwise_plot = stats_popwise_plot

}
