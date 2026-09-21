include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { SPLIT_VCF_BY_CHROM } from "../workflow/utils/SPLIT_VCF_BY_CHROM.nf"
include { RUN_VCF_FILTERING } from "../workflow/run/RUN_VCF_FILTERING.nf"
include { RUN_VCF_FILTERING_POPWISE } from "../workflow/run/RUN_VCF_FILTERING_POPWISE.nf"
include { JOIN_VCF_BY_CHROM } from "../workflow/utils/JOIN_VCF_BY_CHROM.nf"
include { RUN_VCF_QC } from "../workflow/run/RUN_VCF_QC.nf"
include { keyFor } from "../library/filekeys.nf"

def inputChromwise() { file(params.fv_vcf).isDirectory() }
def filterPopwise() { params.fv_popwise as boolean }
def filtersLabel() { keyFor(file(params.fv_flags).name, ["txt"]) }
def publishPath() { "filter_variants/${filtersLabel()}" + (filterPopwise() ? "/popwise" : "") }

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.fv_vcf, params.ref_exclude_coords, genome.chrom_names, true, true)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf, null)

    chroms = inputChromwise()
        ? input.vcf
        : SPLIT_VCF_BY_CHROM(input.vcf_indexed, genome.chrom_names)

    filtered = filterPopwise()
        ? RUN_VCF_FILTERING_POPWISE(chroms, metadata.focal_populations_censuses, params.fv_flags)
        : RUN_VCF_FILTERING(chroms, params.fv_flags)

    concatenated = JOIN_VCF_BY_CHROM(filtered.vcf_filtered, genome.chrom_names, filtersLabel())

    qc = RUN_VCF_QC(
        concatenated.vcf_concat,
        metadata.focal_populations_censuses,
        metadata.population_metadata,
        params.fv_qc_thinning_target,
        params.fv_qc_snpden_binsize,
        genome.chrom_names,
        genome.chrom_labels
    )

    publish:
    filters = filtered.flags
    chrom_vcfs = filtered.vcf_filtered
    concat_vcf = concatenated.vcf_concat
    qc_data = qc.data
    qc_plot = qc.plot
    qc_popwise_data = qc.popwise_data
    qc_popwise_plot = qc.popwise_plot

}

output {

    filters { path "${publishPath()}" }
    concat_vcf { path "${publishPath()}" }
    chrom_vcfs { path "${publishPath()}/chroms" }
    qc_data { path "${publishPath()}/qc/across_pops/data" }
    qc_plot { path "${publishPath()}/qc/across_pops" }
    qc_popwise_data { path "${publishPath()}/qc/within_pops/data" }
    qc_popwise_plot { path "${publishPath()}/qc/within_pops" }

}
