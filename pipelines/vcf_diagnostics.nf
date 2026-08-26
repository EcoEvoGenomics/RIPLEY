include { PARSE_REFERENCE_GENOME } from "../source/nextflow/workflows/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../source/nextflow/workflows/PARSE_VCF.nf"
include { THIN_VCF } from "../source/nextflow/workflows/THIN_VCF.nf"
include { RUN_SNP_DENSITY } from "../source/nextflow/workflows/RUN_SNP_DENSITY.nf"
include { RUN_VCF_STATS } from "../source/nextflow/workflows/RUN_VCF_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.vd_vcf, params.ref_exclude_coords, genome.chrom_names, true)

    RUN_SNP_DENSITY(input.vcf_condensed, params.vd_snpden_binsize, genome.chrom_names, genome.chrom_labels)
    THIN_VCF(input.vcf_annotated, params.vd_thin_to) | RUN_VCF_STATS

    publish:
    snpden_data = RUN_SNP_DENSITY.out.data
    snpden_plot = RUN_SNP_DENSITY.out.plot
    stats_data = RUN_VCF_STATS.out.data
    stats_plot = RUN_VCF_STATS.out.plot

}

output {

    snpden_data { path "vcf_diagnostics" }
    snpden_plot { path "vcf_diagnostics" }
    stats_data { path "vcf_diagnostics" }
    stats_plot { path "vcf_diagnostics" }

}
