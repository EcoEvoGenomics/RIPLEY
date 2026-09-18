include { SAMTOOLS_INDEX as INDEX_CRAM_IN; SAMTOOLS_INDEX as INDEX_CRAM_TMP; SAMTOOLS_INDEX as INDEX_CRAM_OUT } from "../process/samtools.nf"
include { SAMTOOLS_VIEW_TARGETS as SAMTOOLS_PICK_COORDS;  } from "../process/samtools.nf"
include { SAMTOOLS_VIEW_TARGETS as SAMTOOLS_PICK_CHROMS } from "../process/samtools.nf"
include { alphanumericIssue; keyFor } from "../library/filekeys.nf"

workflow PARSE_CRAM {

    take:
    cram_path
    genome_fasta
    genome_fai
    exclude_coords
    chrom_indices
    permit_solo
    permit_dir

    main:
    def permitted_extensions = ["cram"]
    def input_is_solo = (file(cram_path).isFile() && keyFor(file(cram_path).name, permitted_extensions) != null)
    def input_is_dir = file(cram_path).isDirectory()
    if (input_is_solo && input_is_dir) { error("The input path may be interpreted both as file and directory.") }
    if (!(input_is_solo || input_is_dir)) { error("The input path does not exist or is not a directory or CRAM.") }
    if (input_is_solo && !permit_solo) { error("This pipeline cannot process a single CRAM file, only a directory.") }
    if (input_is_dir && !permit_dir) { error("This pipeline cannot process a directory, only a single CRAM file.") }

    cram_ref = genome_fasta.combine(genome_fai)

    if (input_is_dir) {
        cram = Channel.fromPath("${cram_path}/**.cram", checkIfExists: true)
    }

    if (input_is_solo) {
        cram = Channel.fromPath("${cram_path}", checkIfExists: true)
    }

    // The whole name minus its extension is the sample key, and downstream tokenising splits on both underscores and dots
    cram.subscribe { it ->
        def issue = alphanumericIssue(keyFor(it.name, permitted_extensions), "sample key", "Input CRAM ${it.name}")
        if (issue) { error(issue) }
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
            def start = 1
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
