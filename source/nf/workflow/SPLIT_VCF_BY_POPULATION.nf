include { BCFTOOLS_PICK_SAMPLES } from "../process/bcftools.nf"

workflow SPLIT_VCF_BY_POPULATION {

    take:
    vcfs
    population_censuses

    main:
    split_vcfs = population_censuses | combine(vcfs) | BCFTOOLS_PICK_SAMPLES

    emit:
    split_vcfs

}
