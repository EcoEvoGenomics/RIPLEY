include { PLINK_TO_VCF; PLINK_WRITE_SNPLIST; PLINK_EXTRACT_SITES } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

workflow RUN_ADMIXTURE {

    take:
    plinkfiles
    kmin
    kmax
    aim_variance_threshold

    main:
    k_values = Channel.of(kmin..kmax)
    admixture = ADMIXTURE(plinkfiles, k_values)

    snp_list = PLINK_WRITE_SNPLIST(plinkfiles)
    aim_snps = ADMIXTURE_AIMS(admixture.pfile, snp_list, aim_variance_threshold)
    aim_vcfs = PLINK_EXTRACT_SITES(plinkfiles, aim_snps.flatten()) | PLINK_TO_VCF

    emit:
    data = admixture.concat()
    aims = aim_vcfs

}
