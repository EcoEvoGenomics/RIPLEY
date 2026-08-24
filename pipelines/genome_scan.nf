include { PARSE_REFERENCE_GENOME } from "../source/nextflow/workflows/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../source/nextflow/workflows/PARSE_VCF.nf"
include { PARSE_METADATA } from "../source/nextflow/workflows/PARSE_METADATA.nf"
include { RUN_POPGEN_WINDOWS } from "../source/nextflow/workflows/RUN_POPGEN_WINDOWS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix)
    input = PARSE_VCF(params.gs_vcfdir, params.ref_exclude_coords, genome.n_chroms, genome.included, true)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, input.vcf_condensed)

    RUN_POPGEN_WINDOWS(
        input.vcf_annotated,
        metadata.for_samples,
        metadata.focal_populations,
        params.gs_window_size,
        params.gs_step_size,
        params.gs_min_sites
    )

}

output {



}
