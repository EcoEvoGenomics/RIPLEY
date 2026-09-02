include { VCFTOOLS_SNP_DENSITY } from "../process/vcftools.nf"
include { PLOT_VCFTOOLS_SNP_DENSITY } from "../process/plotting.nf"

workflow RUN_SNP_DENSITY {

    take:
    vcf
    bin_size
    chrom_names
    chrom_labels

    main:
    chrom_flag = chrom_names
        .collect()
        .map { chroms -> chroms.join(",") }

    data = VCFTOOLS_SNP_DENSITY(vcf, bin_size)
    plot = PLOT_VCFTOOLS_SNP_DENSITY(data, chrom_flag, chrom_labels)

    emit:
    data = data
    plot = plot
    
}
