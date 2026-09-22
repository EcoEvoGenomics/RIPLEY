include { BEDTOOLS_MAKEWINDOWS } from "../../../common/process/bedtools.nf"
include { BEDTOOLS_INTERSECT_WINDOWS } from "../../process/bedtools.nf"
include { SHAPEIT5_PHASE_COMMON; SHAPEIT5_LIGATE } from "../../process/shapeit5.nf"
include { BCFTOOLS_INDEX } from "../../../common/process/bcftools.nf"
include { BCFTOOLS_BCF_TO_VCF } from "../../process/bcftools.nf"

workflow RUN_VCF_PHASING {

    // Window keys in occupied_windows are <chrom>_<index>, so their leading token matches a chromosome VCF name

    take:
    vcfs_indexed
    genome_index
    window_size
    window_overlap

    main:
    genome_windows = BEDTOOLS_MAKEWINDOWS(genome_index, window_size, window_size - window_overlap).bed_base_zero // Windows must overlap for downstream ligation
    occupied_windows = BEDTOOLS_INTERSECT_WINDOWS(genome_windows.first(), vcfs_indexed.map { vcf, _csi -> vcf }).regions_base_one

    phasing_inputs = occupied_windows
        .splitCsv(sep: "\t")
        .combine(vcfs_indexed)
        .filter { it -> it[0].tokenize("_")[0] == it[2].simpleName }
        .map { it ->
            def name = it[0]
            def region = it[1]
            def vcf = it[2]
            def csi = it[3]
            tuple(name, region, vcf, csi)
        }

    phased_chunks = SHAPEIT5_PHASE_COMMON(phasing_inputs)
        .map { bcf, csi -> tuple(bcf.simpleName.tokenize("_")[0], [bcf, csi]) }
        .groupTuple(by: 0)
        .map { chrom, chunks -> tuple(chrom, chunks.flatten()) }

    phased_vcf = SHAPEIT5_LIGATE(phased_chunks) | BCFTOOLS_BCF_TO_VCF
    phased_vcf_indexed = BCFTOOLS_INDEX(phased_vcf)

    emit:
    vcf = phased_vcf
    vcf_indexed = phased_vcf_indexed

}
