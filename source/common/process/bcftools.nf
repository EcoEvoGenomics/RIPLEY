process BCFTOOLS_INDEX {

    label "BCFTOOLS"

    input:
    path(vcf)

    output:
    tuple path(vcf, includeInputs: true), path("${vcf.name}.csi")

    script:
    """
    bcftools index --threads ${task.cpus} ${vcf}
    """
}

process BCFTOOLS_COUNT_RECORDS {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)

    output:
    path("${vcf.simpleName}.nrecords_chrom.tsv"), emit: per_chrom
    path("${vcf.simpleName}.nrecords.txt"), emit: nrecords

    script:
    """
    bcftools index --stats ${vcf} > ${vcf.simpleName}.nrecords_chrom.tsv
    bcftools index --nrecords ${vcf} > ${vcf.simpleName}.nrecords.txt
    """
}

process BCFTOOLS_PICK_CHROM {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)
    each(chrom)

    output:
    path("${chrom}.vcf.gz")

    script:
    """
    bcftools view \
        --threads ${task.cpus} \
        --regions ${chrom} \
        --output-type z --output ${chrom}_tmp.vcf.gz \
        ${vcf}
    mv ${chrom}_tmp.vcf.gz ${chrom}.vcf.gz
    """
}

process BCFTOOLS_LIST_SAMPLES {

    label "BCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.samples.txt")

    script:
    """
    bcftools query --list-samples ${vcf} > ${vcf.simpleName}.samples.txt
    """
}

process BCFTOOLS_PICK_SAMPLES {

    label "BCFTOOLS"

    input:
    tuple path(sample_list), path(vcf)

    output:
    path("${vcf.simpleName}_${sample_list.simpleName}.vcf.gz"), emit: samples_vcf

    script:
    """
    bcftools view \
        --samples-file ${sample_list} \
        --force-samples \
        --output-type z --output ${vcf.simpleName}_${sample_list.simpleName}.vcf.gz \
        ${vcf}
    """
}

process BCFTOOLS_FILTER_CHROMS {

    label "BCFTOOLS"

    input:
    path(vcf)
    val(keep_chrom_string)

    output:
    path("${vcf.simpleName}.vcf.gz")

    script:
    """
    bcftools view \
        --threads ${task.cpus} \
        --targets ${keep_chrom_string} \
        --output-type z --output ${vcf.simpleName}_tmp.vcf.gz \
        ${vcf}
    mv ${vcf.simpleName}_tmp.vcf.gz ${vcf.simpleName}.vcf.gz
    """
}

process BCFTOOLS_MERGE_VCFS {

    // Inputs must be ordered prior to merging

    label "BCFTOOLS"

    input:
    tuple val(outname), path(vcfs, stageAs: "vcfs/*")

    output:
    path("${outname}.vcf.gz")

    script:
    """
    for vcf in ${vcfs}
    do
        bcftools index --threads ${task.cpus} \${vcf}
        echo "\${vcf}" >> merge.list
    done

    # bcftools merge requires at least two files; single VCF passes through
    if [ \$(wc -l < merge.list) -eq 1 ]
    then
        cp -L \$(cat merge.list) ${outname}.vcf.gz
    else
        bcftools merge \
            --threads ${task.cpus} \
            --file-list merge.list \
            --output-type z --output ${outname}.vcf.gz
    fi
    """
}

process BCFTOOLS_CONCAT_VCFS {

    // Inputs must be ordered prior to concatenation

    label "BCFTOOLS"

    input:
    path(vcfs, stageAs: "vcfs/*")
    val(outname)

    output:
    path("${outname}.vcf.gz")

    script:
    """
    bcftools concat \
        --threads ${task.cpus} \
        --output-type z --output ${outname}.vcf.gz \
        ${vcfs}
    """
}

