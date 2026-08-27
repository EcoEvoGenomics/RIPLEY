include { WRITE_POPULATION_CENSUS } from "../modules/system.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../modules/vcftools.nf"
include { BCFTOOLS_PICK_SAMPLES } from "../modules/bcftools.nf"

workflow SPLIT_VCF {

    take:
    vcfs
    metadata
    populations

    main:
    split_vcfs = WRITE_POPULATION_CENSUS(populations, metadata) | combine(vcfs) | BCFTOOLS_PICK_SAMPLES

    emit:
    split_vcfs

}
