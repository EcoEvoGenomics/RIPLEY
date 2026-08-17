include { PLINK_PCA } from "../modules/plink.nf"

workflow RUN_PCA {

    take:
    bedfiles

    main:
    results = PLINK_PCA(bedfiles)

    emit:
    results

}
