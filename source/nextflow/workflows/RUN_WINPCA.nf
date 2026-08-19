include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../modules/winpca.nf"

workflow RUN_WINPCA {

    take:
    vcf
    genome_index
    metadata

    main:
    repo = GET_WINPCA()
    data = WINPCA_CHROM(repo, metadata, vcf.combine(genome_index.parsed))
    genome_plot = WINPCA_GENOMEPLOT(repo, vcf, metadata, data.data.collect(), genome_index.kept)

    emit:
    data = data.data
    plot = data.plot
        .combine(genome_plot)
    
}
