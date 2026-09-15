include { METADATA_LIST_POPULATION_MEMBERS } from "../process/metadata.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../process/vcftools.nf"
include { BCFTOOLS_PICK_SAMPLES } from "../process/bcftools.nf"

workflow SPLIT_VCF {

    take:
    vcfs
    metadata
    populations

    main:
    split_vcfs = METADATA_LIST_POPULATION_MEMBERS(populations, metadata) | combine(vcfs) | BCFTOOLS_PICK_SAMPLES

    emit:
    split_vcfs

}
