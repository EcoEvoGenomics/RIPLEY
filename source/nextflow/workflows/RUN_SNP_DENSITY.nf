include { VCFTOOLS_SNP_DENSITY; PLOT_VCFTOOLS_SNP_DENSITY } from "../modules/vcftools.nf"

workflow RUN_SNP_DENSITY {

    take:
    vcf
    bin_size
    chr_conversion_table
    plot_chroms

    main:
    data = VCFTOOLS_SNP_DENSITY(vcf, bin_size)
    plot = PLOT_VCFTOOLS_SNP_DENSITY(data, chr_conversion_table, plot_chroms)

    emit:
    data = data
    plot = plot
}
