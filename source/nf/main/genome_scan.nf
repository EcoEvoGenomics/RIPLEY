include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF as PARSE_VCF_SELECTION; PARSE_VCF as PARSE_VCF_STRUCTURE } from "../workflow/PARSE_VCF.nf"
include { PARSE_METADATA as PARSE_METADATA_SELECTION; PARSE_METADATA as PARSE_METADATA_STRUCTURE } from "../workflow/PARSE_METADATA.nf"
include { SPLIT_VCF as SPLIT_VCF_SELECTION; SPLIT_VCF as SPLIT_VCF_STRUCTURE } from "../workflow/SPLIT_VCF.nf"
include { RUN_POPGEN_WINDOWS_SCAN } from "../workflow/RUN_POPGEN_WINDOWS_SCAN.nf"
include { RUN_EHH_SCAN } from "../workflow/RUN_EHH_SCAN.nf"
include { RUN_WINDOWED_PCA_SCAN } from "../workflow/RUN_WINDOWED_PCA_SCAN.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input_selection = PARSE_VCF_SELECTION(params.gs_vcfdir_selection, params.ref_exclude_coords, genome.chrom_names, false, true)
    input_structure = PARSE_VCF_STRUCTURE(params.gs_vcfdir_structure, params.ref_exclude_coords, genome.chrom_names, false, true)
    metadata_selection = PARSE_METADATA_SELECTION(params.metadata, params.focal_populations, input_selection.vcf_condensed)
    metadata_structure = PARSE_METADATA_STRUCTURE(params.metadata, params.focal_populations, input_structure.vcf_condensed)

    population_vcfs_selection = SPLIT_VCF_SELECTION(input_selection.vcf_condensed, metadata_selection.for_samples, metadata_selection.focal_populations)
    population_vcfs_structure = SPLIT_VCF_STRUCTURE(input_structure.vcf_condensed, metadata_structure.for_samples, metadata_structure.focal_populations)

    RUN_POPGEN_WINDOWS_SCAN(
        population_vcfs_selection,
        metadata_selection.for_samples,
        params.gs_window_size,
        params.gs_step_size,
        params.gs_min_sites
    )

    RUN_EHH_SCAN(
        population_vcfs_selection,
        params.gs_window_size,
        params.gs_step_size,
        params.gs_min_sites
    )

    RUN_WINDOWED_PCA_SCAN(
        population_vcfs_structure,
        genome.chrom_indices,
        params.gs_window_size,
        params.gs_step_size
    )

    publish:
    popgen = RUN_POPGEN_WINDOWS_SCAN.out.popgen
    ihs = RUN_EHH_SCAN.out.ihs
    xpehh = RUN_EHH_SCAN.out.xpehh
    wpca = RUN_WINDOWED_PCA_SCAN.out.data

}

output {
    
    popgen { path "genome_scan/popgen" }
    ihs { path "genome_scan/ehh" }
    xpehh { path "genome_scan/ehh" }
    wpca { path "genome_scan/pca" }

}
