include { PLINK_PCA } from "../modules/plink.nf"

workflow RUN_PCA {

    take:
    bedfiles

    main:
    data = PLINK_PCA(bedfiles)

    emit:
    data

}
