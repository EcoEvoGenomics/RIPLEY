include { BCFTOOLS_CONCAT_VCFS } from "../../process/bcftools.nf"
include { tokenCount } from "../../library/filekeys.nf"

workflow JOIN_VCF_BY_CHROM {

    // TAKE -------------------------------------------------
    // vcfs           <- [chr2.vcf.gz, chr1.vcf.gz]
    // ref_chromnames <- [chr1, chr2]
    // new_key        <- variants
    // EMIT -------------------------------------------------
    // vcf_concat     -> variants.vcf.gz
    // ------------------------------------------------------
    //
    // TAKE -------------------------------------------------
    // vcfs           <- [chr2_PopA.vcf.gz, chr1_PopA.vcf.gz,
    //                    chr1_PopB.vcf.gz, chr2_PopB.vcf.gz]
    // ref_chromnames <- [chr1, chr2]
    // new_key        <- variants
    // EMIT -------------------------------------------------
    // vcf_concat     -> [variants_PopA.vcf.gz,
    //                    variants_PopB.vcf.gz]
    // ------------------------------------------------------
    //
    // Always guard with PARSE_VCF upstream to ensure correct
    // order of chromosomes in joined file.
    //
    // Concatenation requires identical samples, so each
    // population is concatenated separately.
    //
    // Inverse of SPLIT_VCF_BY_CHROM.

    take:
    vcfs
    ref_chromnames
    new_key

    main:
    vcfs.map { vcf -> if (tokenCount(vcf.simpleName, "_") > 2) error("Badly tokenized input to JOIN_VCF_BY_CHROM.") }

    ordered = vcfs
        .collect()
        .map { found -> [found] }
        .combine(ref_chromnames.collect().map { names -> [names] })
        .combine(channel.of(new_key))
        .flatMap { i ->
            def found = i[0]
            def chroms = i[1]
            def key = i[2]
            found
                .groupBy { vcf -> vcf.simpleName.tokenize("_")[1] ?: "" }
                .sort { entry -> entry.key }
                .collect { population, group ->
                    tuple(
                        population ? "${key}_${population}" : "${key}",
                        chroms.collect { chrom ->
                            group
                                .findAll { vcf -> vcf.simpleName.tokenize("_")[0] == chrom }
                                .sort { a, b -> a.name <=> b.name }
                        }.flatten()
                    )
                }
        }

    vcf_concat = BCFTOOLS_CONCAT_VCFS(ordered)

    emit:
    vcf_concat = vcf_concat

}
