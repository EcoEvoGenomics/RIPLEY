include { SAMTOOLS_INDEX_CRAM as INDEX_IN; SAMTOOLS_INDEX_CRAM as INDEX_TMP; SAMTOOLS_INDEX_CRAM as INDEX_OUT } from "../process/samtools.nf"
include { SAMTOOLS_EXTRACT_CRAM as SAMTOOLS_PICK_COORDS;  } from "../process/samtools.nf"
include { SAMTOOLS_EXTRACT_CRAM as SAMTOOLS_PICK_CHROMS } from "../process/samtools.nf"

workflow PARSE_CRAM {

    take:
    cram_path
    genome_fasta
    genome_fai
    exclude_coords
    chrom_indices

    main:
    cram = Channel.fromPath("${cram_path}", checkIfExists: true)
    cram_idx = INDEX_IN(cram, genome_fasta, genome_fai)

    exclude_coords ?: null
    if (exclude_coords == null) {
        cram_tmp_idx = cram_idx
    } else {
        cram_tmp = SAMTOOLS_PICK_COORDS(cram_idx, genome_fasta, genome_fai, exclude_coords).drop
        cram_tmp_idx = INDEX_TMP(cram_tmp.keep, genome_fasta, genome_fai)
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
            sort: { idx -> def chrom = idx[0]; chrom }
        )
    
    cram_out = SAMTOOLS_PICK_CHROMS(cram_tmp_idx, genome_fasta, genome_fai, chroms_bed).keep
    cram_out_idx = INDEX_OUT(cram_out, genome_fasta, genome_fai)

    emit:
    parsed = cram_out_idx

}
