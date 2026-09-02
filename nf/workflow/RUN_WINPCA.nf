include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../process/winpca.nf"

workflow RUN_WINPCA {

    take:
    vcf
    genome
    metadata

    main:
    repo = GET_WINPCA()
    data = WINPCA_CHROM(repo, metadata, vcf.combine(genome.index))
    genome_plot = WINPCA_GENOMEPLOT(repo, vcf, metadata, data.data.collect(), genome.included)

    emit:
    data = data.data
    plot = data.plot
        .combine(genome_plot)
    
}
