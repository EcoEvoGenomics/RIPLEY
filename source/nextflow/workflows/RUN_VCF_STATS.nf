include { VCFTOOLS_VCF_STATS } from "../modules/vcftools.nf"
include { PLOT_VCFTOOLS_VCF_STATS } from "../modules/plotting.nf"

workflow RUN_VCF_STATS {
    
    take:
    vcf

    main:
    stats = VCFTOOLS_VCF_STATS(vcf)
    PLOT_VCFTOOLS_VCF_STATS("${launchDir}/source/R/plot_vcf_statistics.R", stats)

    emit:
    plot = PLOT_VCFTOOLS_VCF_STATS.out
    data = stats.collect()

}
