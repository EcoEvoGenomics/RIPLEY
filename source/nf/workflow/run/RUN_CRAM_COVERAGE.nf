include { BEDTOOLS_MAKEWINDOWS } from "../../process/bedtools.nf"
include { SAMTOOLS_BEDCOV } from "../../process/samtools.nf"
include { METADATA_PREPEND_KEY_COLUMN as PREPEND_SAMPLE_COLUMN } from "../../process/metadata.nf"
include { PLOT_SAMTOOLS_BEDCOV } from "../../process/plot.nf"
include { PLOT_SAMTOOLS_BEDCOV as PLOT_SAMTOOLS_BEDCOV_POPWISE } from "../../process/plot.nf"

workflow RUN_CRAM_COVERAGE {

    take:
    cram_indexed
    genome_index
    chrom_names
    chrom_labels
    population_map
    binsize

    main:
    genome_windows = BEDTOOLS_MAKEWINDOWS(genome_index, binsize).bed_base_zero
    bedcov = SAMTOOLS_BEDCOV(cram_indexed.combine(genome_windows))

    bedcov_with_key = bedcov
        .flatten()
        .map { it ->
            def key = it.simpleName
            tuple(key, it)
        } | PREPEND_SAMPLE_COLUMN

    bedcov_combined = bedcov_with_key
        .map { it ->
            def ext = it.extension
            tuple (ext, it)
        }
        .collectFile( { it -> ["coverage.${it[0]}", it[1]]} )

    // Samples absent from the map are non-focal and drop out of the join
    bedcov_popwise = bedcov_with_key
        .map { it -> tuple(it.simpleName, it) }
        .join(population_map)
        .map { _sample, table, population -> tuple(population, table.extension, table) }
        .collectFile( { it -> ["coverage_${it[0]}.${it[1]}", it[2]]} )

    chrom_flag = chrom_names
        .collect()
        .map { chroms -> chroms.join(",") }

    bedcov_plot = PLOT_SAMTOOLS_BEDCOV(bedcov_combined, chrom_flag, chrom_labels, bedcov_combined)
    bedcov_popwise_plot = PLOT_SAMTOOLS_BEDCOV_POPWISE(bedcov_popwise, chrom_flag, chrom_labels, bedcov_combined)

    emit:
    data = bedcov_combined
    plot = bedcov_plot
    popwise_data = bedcov_popwise
    popwise_plot = bedcov_popwise_plot

}
