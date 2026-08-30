include { PLINK_PCA } from "../modules/plink.nf"
include { PLOT_PLINK_PCA } from "../modules/plotting.nf"

workflow RUN_PCA {

    take:
    plinkfiles
    metadata

    main:
    pca = PLINK_PCA(plinkfiles)
    plot = PLOT_PLINK_PCA("${launchDir}/source/R/plot_pca.R", pca.eigenval, pca.eigenvec, metadata)

    emit:
    data = pca.mix()
    plot = plot

}
