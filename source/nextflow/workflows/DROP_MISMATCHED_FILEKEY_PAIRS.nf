workflow DROP_MISMATCHED_FILEKEY_PAIRS {

    // Takes a channel of tuples
    // Removes e.g. [chr1_variantsA.vcf.gz, chr2_variantsB.vcf.gz]
    // Keeps   e.g. [chr1_variantsA.vcf.gz, chr1_variantsB.vcf.gz]
    // Outputs e.g. [chr1, variantsA, variantsB, chr1_variantsA.vcf.gz, chr1_variantsB.vcf.gz]
    
    take:
    file_tuples

    main:
    matched_tuples = file_tuples
        .map { pair ->
            def key_a = file(pair[0]).SimpleName.toString().tokenize("_")[0]
            def key_b = file(pair[1]).SimpleName.toString().tokenize("_")[0]
            def name_a = file(pair[0]).SimpleName.toString().tokenize("_")[1]
            def name_b = file(pair[1]).SimpleName.toString().tokenize("_")[1]
            key_a == key_b
                ? tuple(key_a, name_a, name_b, file(pair[0]), file(pair[1]))
                : null
        }
    
    emit:
    matched_tuples

}
