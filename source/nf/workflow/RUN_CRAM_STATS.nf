include { SAMTOOLS_STATS; SAMTOOLS_COVERAGE; SAMTOOLS_BEDCOV } from "../process/samtools.nf"

workflow RUN_CRAM_STATS {

    take:
    cram_indexed

    main:
    stats = SAMTOOLS_STATS(cram_indexed)
    coverage = SAMTOOLS_COVERAGE(cram_indexed)
    // bedcov = SAMTOOLS_BEDCOV(cram_indexed)

    emit:
    stats = stats
    coverage = coverage
    // bedcov = bedcov

}
