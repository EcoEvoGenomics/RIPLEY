include { VCFTOOLS_CALCULATE_RELATEDNESS } from "../modules/vcftools.nf"
include { PLOT_VCFTOOLS_RELATEDNESS } from "../modules/plotting.nf"

workflow RUN_KINSHIP_ANALYSIS {

    take:
    vcf

    main:
    kinship = VCFTOOLS_CALCULATE_RELATEDNESS(vcf)
    matrix = PLOT_VCFTOOLS_RELATEDNESS(kinship)

    emit:
    data = kinship
    plot = matrix

}
