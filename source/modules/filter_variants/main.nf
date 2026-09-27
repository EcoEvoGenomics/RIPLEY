include { PARSE_REFERENCE_GENOME } from "../../common/workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../../common/workflow/parse/PARSE_VCF.nf"
include { RUN_VCF_FILTERING } from "./workflow/RUN_VCF_FILTERING.nf"
include { JOIN_VCF_BY_CHROM } from "../../common/workflow/utils/JOIN_VCF_BY_CHROM.nf"
include { keyFor } from "../../common/library/filekeys.nf"

def filtersLabel() { keyFor(file(params.flags).name, ["txt"]) }
def publishPath() { "filter_variants/${filtersLabel()}" }

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.vcf, params.ref_exclude_coords, genome.chrom_names)

    filtered = RUN_VCF_FILTERING(input.vcf, params.flags)

    concatenated = JOIN_VCF_BY_CHROM(filtered.vcf_filtered, genome.chrom_names, filtersLabel())

    publish:
    filters = filtered.flags
    chrom_vcfs = filtered.vcf_filtered
    concat_vcf = concatenated.vcf_concat

}

output {

    filters { path "${publishPath()}" }
    concat_vcf { path "${publishPath()}" }
    chrom_vcfs { path "${publishPath()}/chroms" }

}
