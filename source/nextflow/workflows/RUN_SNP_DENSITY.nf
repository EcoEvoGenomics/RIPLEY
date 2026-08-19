include { VCFTOOLS_SNP_DENSITY } from "../modules/vcftools.nf"
include { PLOT_VCFTOOLS_SNP_DENSITY } from "../modules/plotting.nf"

workflow RUN_SNP_DENSITY {

    take:
    vcf
    bin_size
    chrom_conversions
    plot_chroms

    main:
    data = VCFTOOLS_SNP_DENSITY(vcf, bin_size)
    plot = PLOT_VCFTOOLS_SNP_DENSITY("${launchDir}/source/R/plot_snp_density.R", data, plot_chroms, chrom_conversions)

    emit:
    data = data
    plot = plot
    
}
