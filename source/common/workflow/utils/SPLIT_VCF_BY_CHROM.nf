include { BCFTOOLS_PICK_CHROM } from "../../process/bcftools.nf"
include { tokenCount } from "../../library/filekeys.nf"

workflow SPLIT_VCF_BY_CHROM {

    // TAKE ------------------------------------
    // indexed     <- tuple(variants.vcf.gz,
    //                      variants.vcf.gz.csi)
    // chrom_names <- [chr1, chr2]
    // EMIT ------------------------------------
    // split_vcfs  -> [chr1.vcf.gz, chr2.vcf.gz]
    // -----------------------------------------
    //
    // Index used to seek efficiently by region.
    // MUST precede population split as chrom
    // key overwrites the output file name.
    //
    // Inverse of JOIN_VCF_BY_CHROM.

    take:
    indexed
    chrom_names

    main:
    indexed.map { vcf, _index -> if (tokenCount(vcf.simpleName, "_") > 1) error("Badly tokenized input to SPLIT_VCF_BY_CHROM.") }
    split_vcfs = BCFTOOLS_PICK_CHROM(indexed, chrom_names)

    emit:
    split_vcfs

}
