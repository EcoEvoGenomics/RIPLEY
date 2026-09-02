include { BCFTOOLS_INDEX; BCFTOOLS_SAMPLE_VCF } from "../process/bcftools.nf"

workflow THIN_VCF {
    
    take:
    vcf
    n_sites

    main:
    indexed = BCFTOOLS_INDEX(vcf)
    thinned = BCFTOOLS_SAMPLE_VCF(indexed, n_sites)

    emit:
    thinned
    
}
