include { BCFTOOLS_MERGE_VCFS } from "../../process/bcftools.nf"
include { tokenCount } from "../../library/filekeys.nf"

workflow JOIN_VCF_BY_POPULATION {

    // TAKE ------------------------------------
    // vcfs        <- [chr1_PopA.vcf.gz,
    //                 chr1_PopB.vcf.gz,
    //                 chr2_PopA.vcf.gz
    //                 chr2_PopB.vcf.gz]
    // EMIT ------------------------------------
    // joined_vcfs -> [chr1.vcf.gz, chr2.vcf.gz]
    // -----------------------------------------
    //
    // Sample order is not restored, and missing
    // populations are lost.
    //
    // Inverse of SPLIT_VCF_BY_POPULATION.

    take:
    vcfs

    main:
    vcfs.map { vcf -> if (tokenCount(vcf.simpleName, "_") > 2) error("Badly tokenized input to JOIN_VCF_BY_POPULATION.") }

    joined_vcfs = vcfs
        .map { vcf -> tuple(vcf.simpleName.tokenize("_")[0], vcf) }
        .groupTuple(sort: { a, b -> a.name <=> b.name })
        | BCFTOOLS_MERGE_VCFS

    emit:
    joined_vcfs

}
