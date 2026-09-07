include { PLINK_TO_VCF; PLINK_WRITE_SNPLIST; PLINK_EXTRACT_SITES } from "../process/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../process/admixture.nf"
include { BCFTOOLS_VCF_TO_GENOTABLE } from "../process/bcftools.nf"
include { PLOT_ADMIXTURE } from "../process/plotting.nf"

workflow RUN_ADMIXTURE {

    take:
    plinkfiles
    metadata
    kmin
    kmax
    aim_variance_threshold

    main:
    plinkfiles.count().map { n -> n == 1 ?: exit(1, "Can only run ADMIXTURE on one VCF.") }
    
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

    k_min_error = admixture.error
        .reduce { i, j -> j[1] < i[1] ? j : i }
        .map { cv -> def k = cv[0]; k }
    admixture_plot = PLOT_ADMIXTURE(admixture_clusts, k_min_error, metadata)

    aim_snps = ADMIXTURE_AIMS(admixture.alleles, aim_variance_threshold)
    aim_vcfs = PLINK_EXTRACT_SITES(plinkfiles, aim_snps.snpids.flatten()) | PLINK_TO_VCF
    aim_gts = BCFTOOLS_VCF_TO_GENOTABLE(aim_vcfs)
        .map { gt ->
            def stem = gt.tokenize('_')
            def k = stem[-2][1..-1]
            def pops = stem[-1].tokenize('p')
            tuple(k, pops[0], pops[1], gt)
        }

    // TO-DO: Calculate hybrid index and heterozygosity
    // TO-DO: Plot triangle plots

    emit:
    data = admixture.data
    plot = admixture_plot
    clusts = admixture_clusts
    errors = admixture_errors
    aims = aim_gts

}
