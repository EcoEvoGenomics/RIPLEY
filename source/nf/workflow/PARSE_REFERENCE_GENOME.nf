include { asList } from "../library/coerce.nf"
include { alphanumericIssue } from "../library/filekeys.nf"

workflow PARSE_REFERENCE_GENOME {

    take:
    genome_path
    exclude_chroms
    exclude_prefixes
    label_table

    main:
    def exclude = asList(exclude_chroms)
    def prefixes = asList(exclude_prefixes)

    genome = Channel.fromPath("${genome_path}", checkIfExists: true)
    genome_index = Channel.fromPath("${genome_path}.fai", checkIfExists: true)
    genome_annotations = Channel.fromPath("${file(genome_path).parent}/${file(genome_path).baseName}.gff", checkIfExists: true)
    chrom_labels = Channel.fromPath("${label_table}", checkIfExists: true)

    total_chroms = genome_index
        .splitCsv( sep:"\t" )
        .filter { row -> !prefixes.any { prefix -> row[0].toString().startsWith(prefix) } }
        .count()

    index_entries = genome_index
        .splitCsv( sep:"\t" )
        .filter { row -> !prefixes.any { prefix -> row[0].toString().startsWith(prefix) } && !(exclude.contains(row[0])) }
        .ifEmpty { error("There are no contigs in the reference index after filtering.") }

    // Chromosome names survive into filenames downstream, which are tokenised and must be strictly alphanumeric
    index_entries
        .subscribe { row ->
            def issue = alphanumericIssue(row[0], "chromosome", "Reference index ${genome_path}.fai")
            if (issue) { error(issue) }
        }

    chrom_names = index_entries
        .map { row -> row[0] }

    emit:
    fasta = genome
    fai = genome_index
    gff = genome_annotations
    total_chroms = total_chroms
    chrom_indices = index_entries
    chrom_names = chrom_names
    chrom_labels = chrom_labels

}
