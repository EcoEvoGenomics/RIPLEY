include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { SPLIT_VCF_BY_CHROM } from "../workflow/utils/SPLIT_VCF_BY_CHROM.nf"
include { RUN_VCF_FILTERING } from "../workflow/run/RUN_VCF_FILTERING.nf"
include { CONCATENATE_VCFS } from "../workflow/utils/CONCATENATE_VCFS.nf"
include { keyFor } from "../library/filekeys.nf"

nextflow.preview.output = true
def filtersLabel() { keyFor(file(params.fv_filter_flags).name, ["txt"]) }

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.fv_vcf, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, true, true)
    // metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf_condensed, null)

    def input_is_split_by_chrom = file(params.fv_vcf).isDirectory()
    if (input_is_split_by_chrom) {
        chroms = input.vcf_annotated
    } else {
        chroms = SPLIT_VCF_BY_CHROM(input.vcf_annotated_indexed, genome.chrom_names)
    }

    filtered = RUN_VCF_FILTERING(chroms, params.fv_filter_flags)
    concatenated = CONCATENATE_VCFS(filtered.vcf_filtered, genome.chrom_names, Channel.value(filtersLabel()))

    publish:
    filters = filtered.flags
    chrom_vcfs = filtered.vcf_filtered
    concat_vcf = concatenated.vcf

}

output {

    filters { path "filter_variants/${filtersLabel()}" }
    concat_vcf { path "filter_variants/${filtersLabel()}" }
    chrom_vcfs { path "filter_variants/${filtersLabel()}/chroms" }

}
