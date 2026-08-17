include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../modules/winpca.nf"

workflow RUN_WINPCA {

    take:
    vcf
    genome_index
    metadata

    main:
    repo = GET_WINPCA()
    results = WINPCA_CHROM(repo, metadata, vcf.combine(genome_index))
    genome_plot = WINPCA_GENOMEPLOT(repo, vcf, metadata, WINPCA_CHROM.out.data.collect(), genome_index.out.kept)

    emit:
    results = results.out.data
    chrom_plots = results.out.plot
    genome_plot = genome_plot.out
    
}
