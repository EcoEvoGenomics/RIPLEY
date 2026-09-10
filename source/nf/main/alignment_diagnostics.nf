include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/PARSE_METADATA.nf"
include { PARSE_CRAM } from "../workflow/PARSE_CRAM.nf"
include { RUN_CRAM_STATS } from "../workflow/RUN_CRAM_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_CRAM(params.ad_cram, genome.fasta, genome.fai, params.ref_exclude_coords, genome.chrom_indices)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, null, input.parsed)

    stats = RUN_CRAM_STATS(input.parsed, genome.fasta, genome.fai)

    publish:
    data = stats.data

}

output {

    data { path "alignment_diagnostics" }

}
