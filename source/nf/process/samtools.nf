process SAMTOOLS_INDEX {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(ref_fasta), path(ref_fai)

    output:
    tuple \
        path(cram, includeInputs: true), \
        path("${cram}.crai"), \
        path(ref_fasta, includeInputs: true), \
        path(ref_fai, includeInputs: true)

    script:
    """
    samtools index -@ ${task.cpus} ${cram}
    """
}

process SAMTOOLS_VIEW_TARGETS {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai), path(targets)

    output:
    tuple \
        path("keep/${cram.simpleName}.cram"), \
        path(ref_fasta, includeInputs: true), \
        path(ref_fai, includeInputs: true), \
        emit: keep
    tuple \
        path("drop/${cram.simpleName}.cram"), \
        path(ref_fasta, includeInputs: true), \
        path(ref_fai, includeInputs: true), \
        emit: drop

    script:
    """
    mkdir keep
    mkdir drop
    samtools view ${cram} \
        --threads ${task.cpus} \
        --targets-file ${targets} \
        --output keep/${cram.simpleName}.cram \
        --unoutput drop/${cram.simpleName}.cram \
        --reference ${ref_fasta} \
        --cram
    """
}

process SAMTOOLS_STATS {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai)

    output:
    path("${cram.simpleName}.stats")

    script:
    """
    samtools stats ${cram} > ${cram.simpleName}.stats
    """
}

process SAMTOOLS_COVERAGE {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai)

    output:
    path("${cram.simpleName}.coverage")

    script:
    """
    samtools coverage --reference ${ref_fasta} ${cram} > ${cram.simpleName}.coverage
    """
}

process SAMTOOLS_BEDCOV {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai), path(ref_fasta), path(ref_fai), path(windows_bed)

    output:
    path("${cram.simpleName}.bedcov")

    script:
    """
    samtools bedcov ${windows_bed} ${cram} > ${cram.simpleName}.bedcov
    """
}

process SAMTOOLS_FAIDX {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path(fasta, includeInputs: true), path("${fasta.name}.fai")

    script:
    """
    samtools faidx ${fasta}
    """
}

process SAMTOOLS_FAIDX_EXTRACT {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta), path(fai)
    val(region)

    output:
    path("${sample}_${region}.fasta")

    script:
    """
    samtools faidx ${fasta} ${region} > ${sample}_${region}.fasta
    sed -i -e 's/>${region}/>${sample}/g' ${sample}_${region}.fasta
    """
}
