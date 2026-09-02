include { PLINK_TO_VCF; PLINK_WRITE_SNPLIST; PLINK_EXTRACT_SITES } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

workflow RUN_ADMIXTURE {

    take:
    plinkfiles
    metadata
    kmin
    kmax
    aim_variance_threshold

    main:
    k_values = Channel.of(kmin..kmax)

    admixture = ADMIXTURE(plinkfiles, k_values)
    admixture_clusts = admixture.clust.collectFile( name: "admixture.clusts" )
    admixture_errors = admixture.error
        .map { cv ->
            def k = cv[0]
            def cv_error = cv[1]
            "${k}\t${cv_error}\n"
        }
        .collectFile( name: "admixture.errors", sort: { cv -> cv[0]} )

    snp_list = PLINK_WRITE_SNPLIST(plinkfiles)
    aim_snps = ADMIXTURE_AIMS(admixture.pfile, snp_list, aim_variance_threshold)
    aim_vcfs = PLINK_EXTRACT_SITES(plinkfiles, aim_snps.flatten()) | PLINK_TO_VCF

    emit:
    data = admixture.data
    clusts = admixture_clusts
    errors = admixture_errors
    aims = aim_vcfs

}
