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

    cram_ref = genome_fasta.combine(genome_fai)

    cram = Channel.fromPath("${cram_path}/**.cram", checkIfExists: true)

    cram.count().subscribe { n -> if(n < 2) error("The input directory must hold more than one CRAM: ${cram_path}") }

    // The whole name minus its extension is the sample key, and downstream tokenising splits on both underscores and dots
    cram.subscribe { it ->
        def issue = alphanumericIssue(keyFor(it.name, permitted_extensions), "sample key", "Input CRAM ${it.name}")
        if (issue) { error(issue) }
    }

    // Input cram_path is searched recursively, so sample keys can collide across subdirectories.
    cram.map { it -> keyFor(it.name, permitted_extensions) }
        .collect(sort: true)
        .map { keys -> keys.countBy { key -> key }.findAll { _key, n -> n > 1 }.keySet() as List }
        .subscribe { duplicates ->
            if (duplicates) {
                error("The input directory ${cram_path} holds more than one CRAM for sample key or keys: ${duplicates.join(', ')}.")
            }
        }

    cram_idx = cram.combine(cram_ref)
        .map { it ->
            def file = it[0]
            def fasta = it[1]
            def fai = it[2]
            tuple(file, fasta, fai)
        } | INDEX_CRAM_IN

    // combine() rejects a bare path, so the BED must be wrapped to broadcast across every CRAM
    def exclude_bed = exclude_coords ? Channel.value(file(exclude_coords, checkIfExists: true)) : null
    if (exclude_bed == null) {
        cram_tmp = cram_idx
    } else {
        cram_tmp = SAMTOOLS_PICK_COORDS(cram_idx.combine(exclude_bed)).drop | INDEX_CRAM_TMP
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

    cram_out = SAMTOOLS_PICK_CHROMS(cram_tmp.combine(chroms_bed)).keep | INDEX_CRAM_OUT

    emit:
    parsed = cram_out

}
