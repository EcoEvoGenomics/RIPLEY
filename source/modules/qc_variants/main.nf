include { PARSE_REFERENCE_GENOME } from "../../common/workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../../common/workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../../common/workflow/parse/PARSE_VCF.nf"
include { JOIN_VCF_BY_CHROM } from "../../common/workflow/utils/JOIN_VCF_BY_CHROM.nf"
include { RUN_VCF_QC } from "./workflow/RUN_VCF_QC.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.vcf, params.ref_exclude_coords, genome.chrom_names)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf, null)

    // Frequency-based statistics are only meaningful genome-wide
    joined = JOIN_VCF_BY_CHROM(input.vcf, genome.chrom_names, "genome")

    RUN_VCF_QC(
        joined.vcf_concat,
        metadata.focal_populations_censuses,
        metadata.population_metadata,
        params.qc_thinning_target,
        params.qc_snpden_binsize,
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

    data { path "qc_variants/across_pops/data" }
    plot { path "qc_variants/across_pops" }
    popwise_data { path "qc_variants/within_pops/data" }
    popwise_plot { path "qc_variants/within_pops" }

}
