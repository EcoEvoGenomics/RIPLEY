include { PLINK_PCA } from "../process/plink.nf"
include { PLOT_PLINK_PCA } from "../process/plot.nf"

workflow RUN_PCA {

    take:
    plinkfiles
    sample_metadata
    population_metadata
    species_metadata

    main:
    pca = PLINK_PCA(plinkfiles)
    plot = PLOT_PLINK_PCA(pca.eigenval, pca.eigenvec, sample_metadata, population_metadata, species_metadata)

    emit:
    data = pca.mix()
    plot = plot

}
