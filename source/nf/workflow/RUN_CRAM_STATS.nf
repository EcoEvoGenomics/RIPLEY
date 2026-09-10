include { SAMTOOLS_STATS } from "../process/samtools.nf"

workflow RUN_CRAM_STATS {

    take:
    cram
    genome
    index

    main:
    stats = SAMTOOLS_STATS(cram, genome, index)

    emit:
    data = stats.cov.mix(stats.stat)

}
