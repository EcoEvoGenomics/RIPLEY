include { BEDTOOLS_MAKEWINDOWS } from "../process/bedtools.nf"
include { SAMTOOLS_BEDCOV; SAMTOOLS_COVERAGE } from "../process/samtools.nf"
include { METADATA_PREPEND_KEY_COLUMN } from "../process/metadata.nf"

workflow RUN_CRAM_COVERAGE {

    take:
    cram_indexed
    genome_index
    binsize

    main:
    genome_windows = BEDTOOLS_MAKEWINDOWS(genome_index, binsize).bed_base_zero
    bedcov = SAMTOOLS_BEDCOV(cram_indexed.combine(genome_windows))

    bedcov_with_key = bedcov
        .flatten()
        .map { it ->
            def header = "SAMPLE"
            def key = it.simpleName
            tuple(header, key, it)
        } | METADATA_PREPEND_KEY_COLUMN

    bedcov_combined = bedcov_with_key
        .map { it ->
            def ext = it.extension
            tuple (ext, it)
        }
        .collectFile( { it -> ["coverage.${it[0]}", it[1]]},
            skip: 1,
            keepHeader: true
        )

    emit:
    data = bedcov_combined

}
