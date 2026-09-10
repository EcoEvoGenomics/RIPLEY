include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/PARSE_VCF.nf"
include { SPLIT_VCF } from "../workflow/SPLIT_VCF.nf"
include { THIN_VCF; THIN_VCF as THIN_VCF_POPWISE } from "../workflow/THIN_VCF.nf"
include { RUN_SNP_DENSITY } from "../workflow/RUN_SNP_DENSITY.nf"
include { RUN_VCF_STATS; RUN_VCF_STATS as RUN_VCF_STATS_POPWISE } from "../workflow/RUN_VCF_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.vd_vcf, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, true, true)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, input.vcf_condensed)

    RUN_SNP_DENSITY(input.vcf_condensed, params.vd_snpden_binsize, genome.chrom_names, genome.chrom_labels)
    THIN_VCF(input.vcf_annotated, params.vd_thin_to) | RUN_VCF_STATS

    if (metadata.focal_populations_set_by_user) {

        popwise_vcf = SPLIT_VCF(input.vcf_annotated, metadata.for_samples, metadata.focal_populations)
        THIN_VCF_POPWISE(popwise_vcf, params.vd_thin_to) | RUN_VCF_STATS_POPWISE
        
    }

    publish:
    snpden_data = RUN_SNP_DENSITY.out.data
    snpden_plot = RUN_SNP_DENSITY.out.plot
    stats_data = RUN_VCF_STATS.out.data
    stats_plot = RUN_VCF_STATS.out.plot
    popwise_stats_data = RUN_VCF_STATS_POPWISE.out.data
    popwise_stats_plot = RUN_VCF_STATS_POPWISE.out.plot

}

output {

    snpden_data { path "variant_diagnostics" }
    snpden_plot { path "variant_diagnostics" }
    stats_data { path "variant_diagnostics" }
    stats_plot { path "variant_diagnostics" }
    popwise_stats_data { path "variant_diagnostics/pop" }
    popwise_stats_plot { path "variant_diagnostics/pop" }

}
