include { PLINK_PCA } from "../modules/plink.nf"

workflow RUN_PCA {

    take:
    plinkfiles

    main:
    data = PLINK_PCA(plinkfiles)

    emit:
    data

}
