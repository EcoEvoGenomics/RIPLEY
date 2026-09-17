include { BEDTOOLS_MAKEWINDOWS } from "../process/bedtools.nf"
include { SAMTOOLS_BEDCOV; SAMTOOLS_COVERAGE } from "../process/samtools.nf"
include { METADATA_PREPEND_KEY_COLUMN as PREPEND_SAMPLE_COLUMN } from "../process/metadata.nf"
include { PLOT_SAMTOOLS_BEDCOV } from "../process/plot.nf"

workflow RUN_CRAM_COVERAGE {

    take:
    cram_indexed
    genome_index
    chrom_names
    chrom_labels
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
    
    chrom_flag = chrom_names
        .collect()
        .map { chroms -> chroms.join(",") }

    bedcov_plot = PLOT_SAMTOOLS_BEDCOV(bedcov_combined, chrom_flag, chrom_labels)

    emit:
    data = bedcov_combined
    plot = bedcov_plot

}
