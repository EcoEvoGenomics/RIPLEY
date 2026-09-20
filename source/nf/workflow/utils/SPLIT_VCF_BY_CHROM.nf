include { BCFTOOLS_PICK_CHROM } from "../../process/bcftools.nf"

workflow SPLIT_VCF_BY_CHROM {

    take:
    indexed
    chrom_names

    main:
    split_vcfs = BCFTOOLS_PICK_CHROM(indexed, chrom_names)

    emit:
    split_vcfs

}
