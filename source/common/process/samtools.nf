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

