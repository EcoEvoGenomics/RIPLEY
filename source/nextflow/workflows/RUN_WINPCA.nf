include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../modules/winpca.nf"

workflow RUN_WINPCA {

    take:
    vcf
    genome_index
    metadata

    main:
    repo = GET_WINPCA()
    results = WINPCA_CHROM(repo, metadata, vcf.combine(genome_index.parsed))
    genome_plot = WINPCA_GENOMEPLOT(repo, vcf, metadata, results.data.collect(), genome_index.kept)

    emit:
    results = results.data
    chrom_plots = results.plot
    genome_plot = genome_plot
    
}
