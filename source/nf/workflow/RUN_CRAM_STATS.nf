include { SAMTOOLS_STATS; SAMTOOLS_COVERAGE; SAMTOOLS_BEDCOV } from "../process/samtools.nf"
include { BEDTOOLS_MAKEWINDOWS } from "../process/bedtools.nf"

workflow RUN_CRAM_STATS {

    take:
    cram_indexed
    genome_index
    binsize

    main:
    genome_windows = BEDTOOLS_MAKEWINDOWS(genome_index, binsize).bed_base_zero
    
    stats = SAMTOOLS_STATS(cram_indexed)
    coverage_summary = SAMTOOLS_COVERAGE(cram_indexed)
    coverage_details = SAMTOOLS_BEDCOV(cram_indexed.combine(genome_windows))

    emit:
    stats = stats
    coverage_summary = coverage_summary
    coverage_details = coverage_details

}
