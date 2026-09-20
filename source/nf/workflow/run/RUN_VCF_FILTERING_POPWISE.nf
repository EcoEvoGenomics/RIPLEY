include { SPLIT_VCF_BY_POPULATION } from "../utils/SPLIT_VCF_BY_POPULATION.nf"
include { RUN_VCF_FILTERING } from "./RUN_VCF_FILTERING.nf"
include { JOIN_VCF_BY_POPULATION } from "../utils/JOIN_VCF_BY_POPULATION.nf"

workflow RUN_VCF_FILTERING_POPWISE {

    take:
    chrom_vcfs
    population_censuses
    filter_flags_path

    main:
    demerged = SPLIT_VCF_BY_POPULATION(chrom_vcfs, population_censuses)
    filtered = RUN_VCF_FILTERING(demerged, filter_flags_path)
    remerged = JOIN_VCF_BY_POPULATION(filtered.vcf_filtered)

    emit:
    vcf_filtered = remerged
    flags = filtered.flags

}
