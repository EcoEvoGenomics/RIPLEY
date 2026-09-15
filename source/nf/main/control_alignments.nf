include { PARSE_REFERENCE_GENOME } from "../workflow/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../workflow/PARSE_METADATA.nf"
include { PARSE_CRAM } from "../workflow/PARSE_CRAM.nf"
include { RUN_CRAM_STATS } from "../workflow/RUN_CRAM_STATS.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_CRAM(params.ca_cram, genome.fasta, genome.fai, params.ref_exclude_coords, genome.chrom_indices, true, true)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, null, input.parsed)

    stats = RUN_CRAM_STATS(input.parsed)

    if (metadata.focal_populations_set_by_user) {



    }

    publish:
    stats = stats.stats
    coverage = stats.coverage
    // bedcov = stats.bedcov

}

output {

    stats { path "control_alignments/data" }
    coverage { path "control_alignments/data" }
    // bedcov { path "control_alignments/data" }

}
