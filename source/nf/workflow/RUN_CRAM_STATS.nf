include { SAMTOOLS_STAT_CRAM } from "../process/samtools.nf"

workflow RUN_CRAM_STATS {

    take:
    cram_indexed

    main:
    stats = SAMTOOLS_STAT_CRAM(cram_indexed)

    emit:
    data = stats.cov.mix(stats.stat)

}
