include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/PARSE_METADATA.nf"
include { PARSE_CRAM } from "../workflow/PARSE_CRAM.nf"
include { RUN_CRAM_COVERAGE } from "../workflow/RUN_CRAM_COVERAGE.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_CRAM(params.ca_cram, genome.fasta, genome.fai, params.ref_exclude_coords, genome.chrom_indices, true, true)

    coverage = RUN_CRAM_COVERAGE(input.parsed, genome.fai, genome.chrom_names, genome.chrom_labels, params.ca_coverage_binsize)

    publish:
    coverage_data = coverage.data
    coverage_plot = coverage.plot
}

output {

    coverage_data { path "control_alignments/data" }
    coverage_plot { path "control_alignments" }

}
