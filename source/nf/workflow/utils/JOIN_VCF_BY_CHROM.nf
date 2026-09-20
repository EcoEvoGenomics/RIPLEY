include { BCFTOOLS_CONCAT_VCFS } from "../../process/bcftools.nf"

workflow JOIN_VCF_BY_CHROM {

    // SimpleNames or leading tokens of SimpleName are chrom keys.
    // Matching on them here (re)orders by reference genome index.

    take:
    vcfs
    ref_chromnames
    concat_outname

    main:
    ordered = vcfs
        .collect()
        .map { found -> [found] }
        .combine(ref_chromnames.collect().map { names -> [names] })
        .map { i ->
            def found = i[0]
            def chroms = i[1]
            chroms.collect { chrom ->
                found
                    .findAll { vcf -> vcf.simpleName.tokenize("_")[0] == chrom }
                    .sort { a, b -> a.name <=> b.name }
            }.flatten()
        }

    vcf_concat = BCFTOOLS_CONCAT_VCFS(ordered, concat_outname)

    emit:
    vcf = vcf_concat

}
