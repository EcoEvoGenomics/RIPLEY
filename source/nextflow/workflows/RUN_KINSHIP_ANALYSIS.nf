include { VCFTOOLS_CALCULATE_RELATEDNESS; PLOT_VCFTOOLS_RELATEDNESS } from "../modules/vcftools.nf"

workflow RUN_KINSHIP_ANALYSIS {

    take:
    vcf

    main:
    kinship_table = VCFTOOLS_CALCULATE_RELATEDNESS(vcf)
    matrix = PLOT_VCFTOOLS_RELATEDNESS(kinship_table)

    emit:
    table = kinship_table
    matrix = matrix

}
