include { BCFTOOLS_INDEX } from "../../../common/process/bcftools.nf"
include { BCFTOOLS_SAMPLE_VCF } from "../../process/bcftools.nf"

workflow RUN_VCF_THINNING {
    
    take:
    vcf
    n_sites

    main:
    indexed = BCFTOOLS_INDEX(vcf)
    thinned = BCFTOOLS_SAMPLE_VCF(indexed, n_sites)

    emit:
    thinned
    
}
