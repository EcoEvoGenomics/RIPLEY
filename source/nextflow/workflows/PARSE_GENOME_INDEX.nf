workflow PARSE_GENOME_INDEX {

    take:
    index_path
    exclude_chromosomes
    exclude_pattern

    main:
    def pattern = exclude_pattern ?: null
    def exclude = (exclude_chromosomes instanceof List)
        ? exclude_chromosomes as List
        : exclude_chromosomes != null
            ? [exclude_chromosomes]
            : []

    branched_indices = Channel.fromPath(index_path)
        .splitCsv( sep:"\t")
        .map { row -> [row[0], row[1]] }
        .branch { row ->
            excluded: pattern != null && row[0].toString().contains(pattern) || exclude.contains(row[0])
            kept: true
        }

    excluded_index = branched_indices.excluded
    excluded_as_string = excluded_index
        .map { row -> row[0] }
        .collect()
        .map { names -> names.join(",") }

    kept_index = branched_indices.kept
    kept_as_string = kept_index
        .map { row -> row[0] }
        .collect()
        .map { names -> names.join(",") }

    emit:
    parsed = kept_index
    count = kept_index.count()
    excluded = excluded_as_string
    kept = kept_as_string

}
