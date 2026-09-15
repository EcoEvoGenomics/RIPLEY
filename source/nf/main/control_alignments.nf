include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/PARSE_METADATA.nf"
include { PARSE_CRAM } from "../workflow/PARSE_CRAM.nf"
include { RUN_CRAM_COVERAGE } from "../workflow/RUN_CRAM_COVERAGE.nf"
include { COLLATE_MULTIPLE_CRAM_STATS as COLLATE_STATS; COLLATE_MULTIPLE_CRAM_STATS as COLLATE_POPWISE_STATS } from "../workflow/COLLATE_MULTIPLE_CRAM_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_CRAM(params.ca_cram, genome.fasta, genome.fai, params.ref_exclude_coords, genome.chrom_indices, true, true)

    coverage = RUN_CRAM_COVERAGE(input.parsed, genome.fai, params.ca_coverage_binsize)

    publish:
    coverage_data = coverage.data
}

output {

    coverage_data { path "control_alignments/data" }

}
