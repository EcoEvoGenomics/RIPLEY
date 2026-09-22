include { BCFTOOLS_PICK_SAMPLES } from "../../process/bcftools.nf"
include { tokenCount } from "../../library/filekeys.nf"

workflow SPLIT_VCF_BY_POPULATION {

    // TAKE ------------------------------------
    // vcfs       <- [chr1.vcf.gz, chr2.vcf.gz]
    // pop_lists  <- [PopA.list, PopB.list]
    // EMIT ------------------------------------
    // split_vcfs -> [chr1_PopA.vcf.gz,
    //                chr1_PopB.vcf.gz,
    //                chr2_PopA.vcf.gz,
    //                chr2_popB.vcf.gz]
    // -----------------------------------------
    //
    // Inverse of JOIN_VCF_BY_POPULATION.

    take:
    vcfs
    pop_lists

    main:
    vcfs.map { vcf -> if (tokenCount(vcf.simpleName, "_") > 1) error("Badly tokenized input to SPLIT_VCF_BY_POPULATION.") }
    split_vcfs = pop_lists | combine(vcfs) | BCFTOOLS_PICK_SAMPLES

    emit:
    split_vcfs

}
