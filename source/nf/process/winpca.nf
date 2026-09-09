process GET_WINPCA {

    // Downloads Moritz Blumer's WinPCA software
    // See https://github.com/MoritzBlumer/winpca

    label "SYSTEM"

    output:
    path("winpca-1.2.1/*")

    script:
    """
    curl -L https://github.com/MoritzBlumer/winpca/archive/refs/tags/v1.2.1.tar.gz > winpca.tar.gz
    tar -zxvf winpca.tar.gz
    """
}

process WINPCA_PCA_CHROMWISE {

    label "WINPCA"

    // Exits with an unhelpful WinPCA error for VCFs with fewer than 10 000 sites.
    errorStrategy "ignore"

    input:
    path(repo)
    tuple path(vcf), val(chrom), val(chrom_length), val(window_size), val(step_size)

    output:
    path("${vcf.simpleName}*")

    script:
    """
    python3 winpca pca ${vcf.simpleName} ${vcf} ${chrom}:1-${chrom_length} \
        --threads ${task.cpus} \
        --window_size ${window_size} \
        --increment ${step_size} \
        --np
    gunzip ${vcf.simpleName}.pc_1.tsv.gz
    gunzip ${vcf.simpleName}.pc_2.tsv.gz
    gunzip ${vcf.simpleName}.hetp.tsv.gz
    gunzip ${vcf.simpleName}.stat.tsv.gz
    """
}
