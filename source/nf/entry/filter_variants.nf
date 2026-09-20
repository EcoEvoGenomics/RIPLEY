include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { SPLIT_VCF_BY_CHROM } from "../workflow/utils/SPLIT_VCF_BY_CHROM.nf"
include { RUN_VCF_FILTERING } from "../workflow/run/RUN_VCF_FILTERING.nf"
include { JOIN_VCF_BY_CHROM } from "../workflow/utils/JOIN_VCF_BY_CHROM.nf"
include { keyFor } from "../library/filekeys.nf"

def inputChromwise() { file(params.fv_vcf).isDirectory() }
def filtersLabel() { keyFor(file(params.fv_flags).name, ["txt"]) }

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.fv_vcf, params.ref_exclude_coords, genome.total_chroms, genome.chrom_names, true, true)
    // metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf_condensed, null)

    chroms = inputChromwise()
        ? input.vcf_annotated
        : SPLIT_VCF_BY_CHROM(input.vcf_annotated_indexed, genome.chrom_names)

    filtered = RUN_VCF_FILTERING(chroms, params.fv_flags)
    concatenated = JOIN_VCF_BY_CHROM(filtered.vcf_filtered, genome.chrom_names, Channel.value(filtersLabel()))

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
