include { PARSE_REFERENCE_GENOME } from "../source/nextflow/workflows/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../source/nextflow/workflows/PARSE_VCF.nf"
include { PARSE_METADATA } from "../source/nextflow/workflows/PARSE_METADATA.nf"
include { SPLIT_VCF } from "../source/nextflow/workflows/SPLIT_VCF.nf"
include { RUN_POPGEN_WINDOWS } from "../source/nextflow/workflows/RUN_POPGEN_WINDOWS.nf"
include { RUN_EHH_SCAN } from "../source/nextflow/workflows/RUN_EHH_SCAN.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.gs_vcfdir, params.ref_exclude_coords, genome.chrom_names, true)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, input.vcf_condensed)

    population_vcfs = SPLIT_VCF(input.vcf_condensed, metadata.for_samples, metadata.focal_populations)

    RUN_POPGEN_WINDOWS(
        population_vcfs,
        metadata.for_samples,
        params.gs_window_size,
        params.gs_step_size,
        params.gs_min_sites
    )

    RUN_EHH_SCAN(
        population_vcfs,
        params.gs_window_size,
        params.gs_step_size,
        params.gs_min_sites
    )

    publish:
    popgen = RUN_POPGEN_WINDOWS.out
    ihs = RUN_EHH_SCAN.out.ihs
    xpehh = RUN_EHH_SCAN.out.xpehh
}

output {
    
    popgen { path "genome_scan" }
    ihs { path "genome_scan" }
    xpehh { path "genome_scan" }

}
