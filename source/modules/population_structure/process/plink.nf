
process PLINK_INIT_PLINKFILES {

    label "PLINK"

    input:
    path(vcf)
    val(n_chroms)

    output:
    tuple(path("${vcf.simpleName}.bed"), path("${vcf.simpleName}.bim"), path("${vcf.simpleName}.fam"), val(n_chroms))

    script:
    """
    plink \
    --vcf ${vcf} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --double-id \
    --set-missing-var-ids @:# \
    --make-bed --out ${vcf.simpleName}
    """
}

process PLINK_TO_VCF {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    path("${bed.simpleName}.vcf.gz")

    script:
    """
    plink --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --output-chr 'chr26' \
    --recode vcf-iid bgz --out ${bed.simpleName}
    """
}

process PLINK_WRITE_SNPLIST {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    path("${bed.simpleName}.snplist")

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --write-snplist

    mv plink.snplist ${bed.simpleName}.snplist
    """
}

process PLINK_LD_PRUNE {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    val(window_size)
    val(step_size)
    val(r2_threshold)

    output:
    path("pruned.in"), emit: pruned_in
    path("pruned.out"), emit: pruned_out

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --indep-pairwise ${window_size}'kb' ${step_size} ${r2_threshold} \
    --make-bed --out ${bed.simpleName}
    mv ${bed.simpleName}.prune.in pruned.in
    mv ${bed.simpleName}.prune.out pruned.out   
    """
}

process PLINK_EXTRACT_SITES {
    
    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    each(sitelist)

    output:
    tuple \
    path("${bed.simpleName}_${sitelist.simpleName}.bed"), \
    path("${bed.simpleName}_${sitelist.simpleName}.bim"), \
    path("${bed.simpleName}_${sitelist.simpleName}.fam"), \
    val(n_chroms)

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --extract ${sitelist} \
    --make-bed --out ${bed.simpleName}_${sitelist.simpleName}
    """
}

process PLINK_PCA {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    path("${bed.simpleName}.eigenval"), emit: eigenval
    path("${bed.simpleName}.eigenvec"), emit: eigenvec
    path("plink.log"), emit: log

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --pca

    mv plink.eigenval ${bed.simpleName}.eigenval
    awk '{\$2=""; \$1=\$1; print}' plink.eigenvec > ${bed.simpleName}.eigenvec
    """
}
