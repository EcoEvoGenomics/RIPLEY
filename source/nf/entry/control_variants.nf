include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { RUN_VCF_QC } from "../workflow/run/RUN_VCF_QC.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.cv_vcf, params.ref_exclude_coords, genome.chrom_names, true, true)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf_annotated, null)

    RUN_VCF_QC(
        input.vcf_annotated,
        metadata.focal_populations_censuses,
        metadata.population_metadata,
        params.cv_qc_thinning_target,
        params.cv_qc_snpden_binsize,
        genome.chrom_names,
        genome.chrom_labels
    )

    publish:
    data = RUN_VCF_QC.out.data
    plot = RUN_VCF_QC.out.plot
    popwise_data = RUN_VCF_QC.out.popwise_data
    popwise_plot = RUN_VCF_QC.out.popwise_plot

}

output {

    data { path "control_variants/across_pops/data" }
    plot { path "control_variants/across_pops" }
    popwise_data { path "control_variants/within_pops/data" }
    popwise_plot { path "control_variants/within_pops" }

}
