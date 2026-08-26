workflow PARSE_REFERENCE_GENOME {

    take:
    genome_path
    exclude_chroms
    exclude_pattern
    label_table

    main:
    def pattern = exclude_pattern ?: null
    def exclude = (exclude_chroms instanceof List)
        ? exclude_chroms as List
        : exclude_chroms != null
            ? [exclude_chroms]
            : []

    genome = Channel.fromPath("${genome_path}", checkIfExists: true)
    genome_index = Channel.fromPath("${genome_path}.fai", checkIfExists: true)
    genome_annotations = Channel.fromPath("${file(genome_path).parent}/${file(genome_path).baseName}.gff", checkIfExists: true)
    chrom_labels = Channel.fromPath("${label_table}", checkIfExists: true)

    index_entries = genome_index
        .splitCsv( sep:"\t" )
        .filter { row -> !(pattern != null && row[0].toString().contains(pattern)) && !(exclude.contains(row[0])) }
        .ifEmpty { exit(1, "There are no contigs in the reference index after filtering.") }

    chrom_names = index_entries
        .map { row -> row[0] }
        
    emit:
    fasta = genome
    fai = genome_index
    gff = genome_annotations
    chrom_indices = index_entries
    chrom_names = chrom_names
    chrom_labels = chrom_labels

}
