include { BCFTOOLS_CONCAT_VCFS } from "../../process/bcftools.nf"

workflow CONCATENATE_VCFS {

    take:
    vcfs
    ref_chromnames
    concat_outname

    main:
    // VCF simpleNames are chrom keys, so exact match here orders by reference .fai
    ordered = vcfs
        .collect()
        .map { found -> [found] }
        .combine(ref_chromnames.collect().map { names -> [names] })
        .map { i ->
            def found = i[0]
            def chroms = i[1]
            chroms.collect { chrom -> found.find { vcf -> vcf.simpleName == chrom } }.findAll()
        }

    concatenated_vcf = BCFTOOLS_CONCAT_VCFS(ordered, concat_outname)

    emit:
    vcf = concatenated_vcf

}
