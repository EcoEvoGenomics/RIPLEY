include { PARSE_REFERENCE_GENOME } from "../../common/workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../../common/workflow/parse/PARSE_VCF.nf"
include { RUN_VCF_FILTERING } from "./workflow/RUN_VCF_FILTERING.nf"
include { RUN_VCF_FILTERING_VARIANT_INVARIANT } from "./workflow/RUN_VCF_FILTERING_VARIANT_INVARIANT.nf"
include { keyFor } from "../../common/library/filekeys.nf"

def separateInvariants() { params.flagsInvariants != null }
def filtersLabel() { keyFor(file(params.flagsVariants).name, ["txt"]) }
def publishPath() { "filter_variants/${filtersLabel()}" }

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.vcf, params.ref_exclude_coords, genome.chrom_names)

    filtered = separateInvariants()
        ? RUN_VCF_FILTERING_VARIANT_INVARIANT(input.vcf, params.flagsVariants, params.flagsInvariants)
        : RUN_VCF_FILTERING(input.vcf, params.flagsVariants)

    publish:
    filters = filtered.flags
    chrom_vcfs = filtered.vcf_filtered

}

output {

    filters { path "${publishPath()}" }
    chrom_vcfs { path "${publishPath()}/chroms" }

}
