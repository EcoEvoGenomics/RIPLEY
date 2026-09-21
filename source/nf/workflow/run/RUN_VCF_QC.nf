include { SPLIT_VCF_BY_POPULATION } from "../utils/SPLIT_VCF_BY_POPULATION.nf"
include { RUN_VCF_THINNING; RUN_VCF_THINNING as RUN_VCF_THINNING_POPWISE } from "./RUN_VCF_THINNING.nf"
include { RUN_SNP_DENSITY; RUN_SNP_DENSITY as RUN_SNP_DENSITY_POPWISE } from "./RUN_SNP_DENSITY.nf"
include { RUN_VCF_STATS; RUN_VCF_STATS as RUN_VCF_STATS_POPWISE } from "./RUN_VCF_STATS.nf"
include { COLLATE_POPWISE_VCF_STATS } from "../utils/COLLATE_POPWISE_VCF_STATS.nf"

workflow RUN_VCF_QC {

    take:
    vcf
    population_censuses
    population_metadata
    thin_to
    snpden_binsize
    chrom_names
    chrom_labels

    main:
    snpden = RUN_SNP_DENSITY(vcf, snpden_binsize, chrom_names, chrom_labels)
    stats = RUN_VCF_THINNING(vcf, thin_to) | RUN_VCF_STATS

    popwise_vcf = SPLIT_VCF_BY_POPULATION(vcf, population_censuses)

    popwise_snpden = RUN_SNP_DENSITY_POPWISE(popwise_vcf, snpden_binsize, chrom_names, chrom_labels)
    popwise_stats = RUN_VCF_THINNING_POPWISE(popwise_vcf, thin_to) | RUN_VCF_STATS_POPWISE
    popwise_collated = COLLATE_POPWISE_VCF_STATS(popwise_stats.data, population_metadata)

    emit:
    data = stats.data.mix(snpden.data)
    plot = stats.plot.mix(snpden.plot)
    popwise_data = popwise_collated.data.mix(popwise_snpden.data)
    popwise_plot = popwise_collated.plot.mix(popwise_snpden.plot)

}
