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
    input_selection = PARSE_VCF_SELECTION(params.sl_vcfdir_selection, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, false, true)
    input_structure = PARSE_VCF_STRUCTURE(params.sl_vcfdir_structure, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, false, true)
    metadata_selection = PARSE_METADATA_SELECTION(params.metadata, params.focal_populations, input_selection.vcf_condensed, null)
    metadata_structure = PARSE_METADATA_STRUCTURE(params.metadata, params.focal_populations, input_structure.vcf_condensed, null)

    popwise_vcf_selection = SPLIT_VCF_SELECTION(input_selection.vcf_condensed, metadata_selection.for_samples, metadata_selection.focal_populations)
    popwise_vcf_structure = SPLIT_VCF_STRUCTURE(input_structure.vcf_condensed, metadata_structure.for_samples, metadata_structure.focal_populations)

    RUN_POPGEN_WINDOWS_SCAN(
        popwise_vcf_selection,
        metadata_selection.for_samples,
        params.sl_window_size,
        params.sl_step_size,
        params.sl_min_sites
    )

    RUN_EHH_SCAN(
        popwise_vcf_selection,
        params.sl_window_size,
        params.sl_step_size,
        params.sl_min_sites
    )

    RUN_WINDOWED_PCA_SCAN(
        popwise_vcf_structure,
        genome.chrom_indices,
        params.sl_window_size,
        params.sl_step_size
    )

    publish:
    popgen = RUN_POPGEN_WINDOWS_SCAN.out.popgen
    ihs = RUN_EHH_SCAN.out.ihs
    xpehh = RUN_EHH_SCAN.out.xpehh
    wpca = RUN_WINDOWED_PCA_SCAN.out.data

}

output {
    
    popgen { path "selection_scan/popgen" }
    ihs { path "selection_scan/ehh" }
    xpehh { path "selection_scan/ehh" }
    wpca { path "selection_scan/pca" }

}
