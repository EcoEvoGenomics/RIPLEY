include { SAMTOOLS_INDEX as INDEX_CRAM_IN; SAMTOOLS_INDEX as INDEX_CRAM_TMP; SAMTOOLS_INDEX as INDEX_CRAM_OUT } from "../../process/samtools.nf"
include { SAMTOOLS_VIEW_TARGETS as SAMTOOLS_PICK_COORDS;  } from "../../process/samtools.nf"
include { SAMTOOLS_VIEW_TARGETS as SAMTOOLS_PICK_CHROMS } from "../../process/samtools.nf"
include { alphanumericIssue; keyFor } from "../../library/filekeys.nf"

workflow PARSE_CRAM {

    take:
    cram_path
    genome_fasta
    genome_fai
    exclude_coords
    chrom_indices

    main:
    def permitted_extensions = ["cram"]

    def caller_input = file(cram_path)
    if (!(caller_input instanceof Path)) { error("The input path must name a directory, not a glob pattern: ${cram_path}") }
    def input_is_dir = caller_input.isDirectory()
    if (!input_is_dir) { error("The input path does not exist or is not a directory: ${cram_path}") }

    crams = Channel.fromPath("${cram_path}/**.cram", checkIfExists: true)

    crams.count().subscribe { n -> if(n < 2) error("The input directory must hold more than one CRAM: ${cram_path}") }
    crams.subscribe { it ->
        // Filnames minus extensions are sample keys. Downstream tokenising relies on strictly alphanumeric keys.
        def issue = alphanumericIssue(keyFor(it.name, permitted_extensions), "sample key", "Input CRAM ${it.name}")
        if (issue) { error(issue) }
    }

    // Input cram_path is searched recursively, so sample keys can collide across subdirectories.
    crams.map { it -> keyFor(it.name, permitted_extensions) }
        .collect(sort: true)
        .map { keys -> keys.countBy { key -> key }.findAll { _key, n -> n > 1 }.keySet() as List }
        .subscribe { duplicates ->
            if (duplicates) {
                error("The input directory ${cram_path} holds more than one CRAM for sample key or keys: ${duplicates.join(', ')}.")
            }
        }

    crams_indexed = crams.combine(genome_fasta.combine(genome_fai))
        .map { it ->
            def file = it[0]
            def fasta = it[1]
            def fai = it[2]
            tuple(file, fasta, fai)
        } | INDEX_CRAM_IN

    def exclude_bed = exclude_coords ? Channel.value(file(exclude_coords, checkIfExists: true)) : null
    if (exclude_bed == null) {
        crams_tmp = crams_indexed
    } else {
        crams_tmp = SAMTOOLS_PICK_COORDS(crams_indexed.combine(exclude_bed)).drop | INDEX_CRAM_TMP
    }

    chroms_bed = chrom_indices
        .map { idx ->
            def chrom = idx[0]
            def start = 0 // Bed file is 0-based
            def end = idx[1]
            "${chrom}\t${start}\t${end}\n"
        }
        .collectFile(
            name: "chroms.bed",
            sort: { line -> line.tokenize("\t")[0] } // Entries are BED lines
        )

    crams_parsed = SAMTOOLS_PICK_CHROMS(crams_tmp.combine(chroms_bed)).keep | INDEX_CRAM_OUT

    emit:
    parsed = crams_parsed

}
