include { PARSE_GENOME_INDEX } from "../source/nextflow/workflows/PARSE_GENOME_INDEX.nf"
include { PARSE_VCF_FILE } from "../source/nextflow/workflows/PARSE_VCF_FILE.nf"
include { THIN_VCF } from "../source/nextflow/workflows/THIN_VCF.nf"
include { RUN_SNP_DENSITY } from "../source/nextflow/workflows/RUN_SNP_DENSITY.nf"
include { RUN_VCF_STATS } from "../source/nextflow/workflows/RUN_VCF_STATS.nf"


nextflow.preview.output = true

workflow {
    main:
    index = PARSE_GENOME_INDEX(params.ref_genome_index, params.ref_exclude_chroms, params.ref_exclude_prefix)
    input = PARSE_VCF_FILE(params.vd_vcf, index.count, index.excluded)
    vcf_annot = input.vcf_annot
    vcf = input.vcf
    
    downsampled = THIN_VCF(vcf_annot, params.vd_thin_to)

    RUN_SNP_DENSITY(vcf, params.vd_snpden_binsize, params.ref_chroms_renamed, index.kept)
    RUN_VCF_STATS(downsampled)

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
