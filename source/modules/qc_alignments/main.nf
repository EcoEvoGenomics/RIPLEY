include { PARSE_REFERENCE_GENOME } from "../../common/workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_METADATA } from "../../common/workflow/parse/PARSE_METADATA.nf"
include { PARSE_CRAM } from "../../common/workflow/parse/PARSE_CRAM.nf"
include { RUN_CRAM_STATS } from "./workflow/RUN_CRAM_STATS.nf"
include { RUN_CRAM_COVERAGE } from "./workflow/RUN_CRAM_COVERAGE.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_CRAM(params.cram, genome.fasta, genome.fai, params.ref_exclude_coords, genome.chrom_indices, true, true)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, null, input.parsed.map { cram_set -> cram_set[0] })

    stats = RUN_CRAM_STATS(input.parsed, metadata.focal_population_map, metadata.population_metadata)
    coverage = RUN_CRAM_COVERAGE(input.parsed, genome.fai, genome.chrom_names, genome.chrom_labels, metadata.focal_population_map, params.coverage_binsize)

    publish:
    stats_data = stats.data
    stats_plot = stats.plot
    stats_popwise_data = stats.popwise_data
    stats_popwise_plot = stats.popwise_plot
    coverage_data = coverage.data
    coverage_plot = coverage.plot
    coverage_popwise_data = coverage.popwise_data
    coverage_popwise_plot = coverage.popwise_plot
}

output {

    stats_data { path "qc_alignments/data" }
    stats_plot { path "qc_alignments" }
    stats_popwise_data { path "qc_alignments/data" }
    stats_popwise_plot { path "qc_alignments" }
    coverage_data { path "qc_alignments/data" }
    coverage_plot { path "qc_alignments" }
    coverage_popwise_data { path "qc_alignments/data" }
    coverage_popwise_plot { path "qc_alignments" }

}
