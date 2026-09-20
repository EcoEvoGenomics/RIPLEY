include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { SPLIT_VCF_BY_POPULATION } from "../workflow/utils/SPLIT_VCF_BY_POPULATION.nf"
include { THIN_VCF; THIN_VCF as THIN_VCF_POPWISE } from "../workflow/run/THIN_VCF.nf"
include { RUN_SNP_DENSITY } from "../workflow/run/RUN_SNP_DENSITY.nf"
include { RUN_VCF_STATS; RUN_VCF_STATS as RUN_VCF_STATS_POPWISE } from "../workflow/run/RUN_VCF_STATS.nf"
include { COLLATE_POPWISE_VCF_STATS } from "../workflow/utils/COLLATE_POPWISE_VCF_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.cv_vcf, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, true, true)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf_condensed, null)

    RUN_SNP_DENSITY(input.vcf_condensed, params.cv_snpden_binsize, genome.chrom_names, genome.chrom_labels)
    THIN_VCF(input.vcf_annotated, params.cv_thin_to) | RUN_VCF_STATS
    
    popwise_vcf = SPLIT_VCF_BY_POPULATION(input.vcf_annotated, metadata.focal_populations_censuses)
    popwise_stats = THIN_VCF_POPWISE(popwise_vcf, params.cv_thin_to) | RUN_VCF_STATS_POPWISE
    COLLATE_POPWISE_VCF_STATS(popwise_stats.data, metadata.population_metadata)

    publish:
    snpden_data = RUN_SNP_DENSITY.out.data
    snpden_plot = RUN_SNP_DENSITY.out.plot
    stats_data = RUN_VCF_STATS.out.data
    stats_plot = RUN_VCF_STATS.out.plot
    popwise_stats_data = COLLATE_POPWISE_VCF_STATS.out.data
    popwise_stats_plot = COLLATE_POPWISE_VCF_STATS.out.plot

}

output {

    snpden_data { path "control_variants/data" }
    stats_data { path "control_variants/data" }
    popwise_stats_data { path "control_variants/data" }
    snpden_plot { path "control_variants" }
    stats_plot { path "control_variants" }
    popwise_stats_plot { path "control_variants" }

}
