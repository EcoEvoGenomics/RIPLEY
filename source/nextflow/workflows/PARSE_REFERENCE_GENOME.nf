workflow PARSE_REFERENCE_GENOME {

    take:
    genome_path
    exclude_chroms
    exclude_pattern

    main:
    def pattern = exclude_pattern ?: null
    def exclude = (exclude_chroms instanceof List)
        ? exclude_chroms as List
        : exclude_chroms != null
            ? [exclude_chroms]
            : []

    genome = Channel.fromPath("${genome_path}", checkIfExists: true)
    genome_annotations = Channel.fromPath("${genome_path.parent}/${genome_path.simpleName}.gff", checkIfExists: true)

    genome_index = Channel.fromPath("${genome_path}.fai", checkIfExists: true)
        .splitCsv( sep:"\t")
        .map { row -> [row[0], row[1]] }
        .branch { row ->
            excluded: pattern != null && row[0].toString().contains(pattern) || exclude.contains(row[0])
            included: true
        }

    excluded_index = genome_index.excluded
    included_index = genome_index.included
    
    included_index
        .ifEmpty { exit(1, "There are no contigs in the reference index after filtering.") }

    excluded_as_string = excluded_index
        .map { row -> row[0] }
        .collect()
        .map { names -> names.join(",") }

    included_as_string = included_index
        .map { row -> row[0] }
        .collect()
        .map { names -> names.join(",") }

    emit:
    genome = genome
    annotations = genome_annotations
    index = included_index
    n_chroms = included_index.count()
    excluded = excluded_as_string
    included = included_as_string

}
