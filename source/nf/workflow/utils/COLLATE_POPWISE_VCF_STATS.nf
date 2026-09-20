include { METADATA_PREPEND_KEY_COLUMN_WITH_HEADER as PREPEND_POP_COLUMN } from "../../process/metadata.nf"
include { PLOT_VCFTOOLS_VCF_STATS_POPWISE } from "../../process/plot.nf"

workflow COLLATE_POPWISE_VCF_STATS {

    take:
    vcf_stats
    population_metadata

    main:
    with_key = vcf_stats
        .flatten()
        .map { it ->
            def population = it.simpleName.tokenize("_")[1]
            tuple("POP", population, it)
        } \
        | PREPEND_POP_COLUMN
    
    // Stem here is original VCF, either single input or one of directory
    across_stems = with_key
        .flatten()
        .map { it ->
            def ext = it.extension
            def stem = it.simpleName.tokenize("_")[0]
            tuple(stem, ext, it)
        }
        .collectFile( { it -> ["${it[0]}_popwise.${it[1]}", it[2]] },
            skip: 1,
            keepHeader: true
        )

    // The plot script reads one file per statistic, so each stem must be plotted separately rather than all stems at once
    per_stem = across_stems
        .map { it -> tuple(it.simpleName.tokenize("_")[0], it) }
        .groupTuple()
        .map { it -> it[1] }

    plot = PLOT_VCFTOOLS_VCF_STATS_POPWISE(per_stem, population_metadata.first())

    emit:
    data = across_stems.collect()
    plot = plot

}
