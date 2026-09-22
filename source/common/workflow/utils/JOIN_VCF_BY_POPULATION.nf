include { BCFTOOLS_MERGE_VCFS } from "../../process/bcftools.nf"

workflow JOIN_VCF_BY_POPULATION {

    // Populations are always added as suffixes, so joining
    // on the leading key restores the original file name

    take:
    vcfs

    main:
    joined_vcfs = vcfs
        .map { vcf -> tuple(vcf.simpleName.tokenize("_")[0], vcf) }
        .groupTuple(sort: { a, b -> a.name <=> b.name })
        | BCFTOOLS_MERGE_VCFS

    emit:
    joined_vcfs

}
