include { VCFTOOLS_VCF_STATS } from "../modules/vcftools.nf"
include { PLOT_VCFTOOLS_VCF_STATS } from "../modules/plotting.nf"

workflow RUN_VCF_STATS {
    
    take:
    vcf

    main:
    VCFTOOLS_VCF_STATS(vcf)
    frq = VCFTOOLS_VCF_STATS.out.frq
    idepth = VCFTOOLS_VCF_STATS.out.idepth
    imiss = VCFTOOLS_VCF_STATS.out.imiss
    ldepth_mean = VCFTOOLS_VCF_STATS.out.ldepth_mean
    lqual = VCFTOOLS_VCF_STATS.out.lqual
    lmiss = VCFTOOLS_VCF_STATS.out.lmiss
    het = VCFTOOLS_VCF_STATS.out.het
    hwe = VCFTOOLS_VCF_STATS.out.hwe
    PLOT_VCFTOOLS_VCF_STATS(vcf, frq, idepth, imiss, ldepth_mean, lqual, lmiss, het, hwe)

    emit:
    plot = PLOT_VCFTOOLS_VCF_STATS.out
    data = frq
        .combine(idepth)
        .combine(imiss)
        .combine(ldepth_mean)
        .combine(lqual)
        .combine(lmiss)
        .combine(het)
        .combine(hwe)

}
