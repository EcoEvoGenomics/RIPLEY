include { PLINK_TO_VCF; PLINK_WRITE_SNPLIST; PLINK_EXTRACT_SITES } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

workflow RUN_ADMIXTURE {

    take:
    bedfiles
    kmin
    kmax
    aim_variance_threshold

    main:
    k_values = Channel.of(kmin..kmax)
    admix = ADMIXTURE(bedfiles, k_values)

    snp_list = PLINK_WRITE_SNPLIST(bedfiles)
    aim_snps = ADMIXTURE_AIMS(admix.pfile, snp_list, aim_variance_threshold)
    vcfs = PLINK_EXTRACT_SITES(bedfiles, aim_snps.flatten()) | PLINK_TO_VCF

    emit:
    data = admix.concat()
    aims = vcfs

}
