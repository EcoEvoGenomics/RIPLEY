include { VCFTOOLS_CALCULATE_RELATEDNESS } from "../process/vcftools.nf"
include { PLOT_VCFTOOLS_RELATEDNESS } from "../process/plot.nf"

workflow RUN_KINSHIP_ANALYSIS {

    take:
    vcf
    sample_metadata
    population_metadata
    species_metadata

    main:
    kinship = VCFTOOLS_CALCULATE_RELATEDNESS(vcf)
    matrix = PLOT_VCFTOOLS_RELATEDNESS(
        file("${moduleDir}/../library/plot_kinship.R", checkIfExists: true),
        kinship, sample_metadata, population_metadata, species_metadata
    )

    emit:
    data = kinship
    plot = matrix

}
