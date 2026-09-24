include { asList } from "../../library/coerce.nf"
include { alphanumericIssue } from "../../library/filekeys.nf"

workflow PARSE_REFERENCE_GENOME {

    take:
    genome_path
    ploidy_path
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
    ploidy = (ploidy_path != null) ? Channel.fromPath("${ploidy_path}", checkIfExists: true) : Channel.empty()
    
    // Assumption: unprefixed contigs = chromosomes; prefixed contigs = unattached scaffolds
    unprefixed_index_entries = genome_index
        .splitCsv( sep:"\t" )
        .filter { row -> !prefixes.any { prefix -> row[0].toString().startsWith(prefix) } }

    // Chromosome names survive into filenames downstream, which are tokenised and must be strictly alphanumeric
    unprefixed_index_entries
        .subscribe { row ->
            def issue = alphanumericIssue(row[0], "chromosome", "Reference index ${genome_path}.fai")
            if (issue) { error(issue) }
        }

    unprefixed_contigs = unprefixed_index_entries
        .map { row -> row[0] }
        .collect(sort: true)

    postfilter_index_entries = unprefixed_index_entries
        .filter { row -> !(exclude.contains(row[0])) }
        .ifEmpty { error("There are no contigs in the reference index after filtering.") }

    ploidy_entries = ploidy
        .splitCsv( sep:"\t" )

    // Ignore wildcards, and accept a ploidy file may legitimately cover excluded chromosomes (but not unattached scaffolds)
    ploidy_entries
        .filter { row -> row[0] != "*" }
        .map { row -> row[0] }
        .combine(unprefixed_contigs.map { names -> [names.unique()] })
        .filter { i -> i[0] !in i[1] }
        .map { i -> i[0] }
        .collect(sort: true)
        .subscribe { absent ->
            if (absent) {
                error("Ploidy file ${ploidy_path} names contig(s) absent from ${genome_path}.fai: ${absent.unique().join(', ')}.")
            }
        }

    ploidy_sexes = ploidy_entries
        .map { row -> row[3] }
        .collect(sort: true)
        .map { sexes -> [sexes.unique()] }

    emit:
    fasta = genome
    fai = genome_index
    gff = genome_annotations
    total_chroms = unprefixed_contigs.map { names -> names.size() }
    chrom_indices = postfilter_index_entries
    chrom_names = postfilter_index_entries.map { row -> row[0] }
    chrom_labels = chrom_labels
    ploidy = ploidy
    ploidy_sexes = ploidy_sexes

}
