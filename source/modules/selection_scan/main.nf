include { PARSE_REFERENCE_GENOME } from "../../common/workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../../common/workflow/parse/PARSE_VCF.nf"
include { PARSE_METADATA } from "../../common/workflow/parse/PARSE_METADATA.nf"
include { SPLIT_VCF_BY_POPULATION } from "../../common/workflow/utils/SPLIT_VCF_BY_POPULATION.nf"
include { RUN_VCF_PHASING } from "./workflow/RUN_VCF_PHASING.nf"
include { RUN_PIXY_SCANS } from "./workflow/RUN_PIXY_SCANS.nf"
include { RUN_EHH_SCAN } from "./workflow/RUN_EHH_SCAN.nf"
include { RUN_DCMS } from "./workflow/RUN_DCMS.nf"

def verifyScanGrid(step_size, window_size) {
    def step_is_immobile = (step_size as int) < 1
    def step_divides_window = (window_size as int) % (step_size as int) == 0
    if (step_is_immobile || !step_divides_window) {
        error("scan_step_size (${step_size}) must be a positive integer factor of scan_window_size (${window_size}).")
    }
}

workflow {

    main:
    verifyScanGrid(params.scan_step_size, params.scan_window_size)
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)

    input = PARSE_VCF(params.vcf, params.ref_exclude_coords, genome.chrom_names)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf, null)
    phased_vcf = RUN_VCF_PHASING(input.vcf_indexed, genome.fai, params.phase_window_size, params.phase_window_overlap)
    popwise_vcf = SPLIT_VCF_BY_POPULATION(phased_vcf.vcf, metadata.focal_populations_censuses)

    RUN_PIXY_SCANS(
        input.vcf_indexed,
        metadata.focal_population_map,
        genome.fai,
        params.pixy_stats,
        params.scan_window_size,
        params.scan_step_size
    )

    RUN_EHH_SCAN(
        popwise_vcf,
        params.scan_window_size,
        params.scan_step_size,
        params.scan_min_sites
    )

    RUN_DCMS(
        RUN_PIXY_SCANS.out.dcms_formatted.mix(RUN_EHH_SCAN.out.dcms_formatted),
        metadata.population_metadata,
        genome.chrom_labels
    )

    publish:
    pixy = RUN_PIXY_SCANS.out.pixy
    ihs = RUN_EHH_SCAN.out.ihs
    xpehh = RUN_EHH_SCAN.out.xpehh
    dcms = RUN_DCMS.out.data
    dcms_plot_main = RUN_DCMS.out.plot

}

output {

    pixy { path "selection_scan/pixy" }
    ihs { path "selection_scan/ehh" }
    xpehh { path "selection_scan/ehh" }
    dcms { path "selection_scan/dcms" }
    dcms_plot_main { path "selection_scan/dcms" }

}
