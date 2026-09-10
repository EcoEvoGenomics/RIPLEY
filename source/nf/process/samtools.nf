process SAMTOOLS_INDEX_CRAM {

    label "SAMTOOLS"

    input:
    path(cram)
    path(ref_fasta)
    path(ref_fai)

    output:
    tuple path(cram, includeInputs: true), path("${cram}.crai")

    script:
    """
    samtools index -@ ${task.cpus} ${cram}
    """
}

process SAMTOOLS_EXTRACT_CRAM {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai)
    path(ref_fasta)
    path(ref_fai)
    path(targets)

    output:
    path("keep/${cram.simpleName}.cram"), emit: keep
    path("drop/${cram.simpleName}.cram"), emit: drop

    script:
    """
    mkdir keep
    mkdir drop
    samtools view ${cram} \
        --threads ${task.cpus} \
        --targets-file ${targets} \
        --unoutput drop/${cram.simpleName}.cram \
        --output keep/${cram.simpleName}.cram \
        --reference ${ref_fasta} \
        --cram
    """
}

process SAMTOOLS_INDEX_FASTA {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path(fasta, includeInputs: true), path("${fasta.name}.fai"), emit: indexed_fasta

    script:
    """
    samtools faidx ${fasta}
    """
}

process SAMTOOLS_EXTRACT_FASTA {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta), path(fai)
    val(region)

    output:
    path("${sample}_${region}.fasta"), emit: extracted

    script:
    """
    samtools faidx ${fasta} ${region} > ${sample}_${region}.fasta
    sed -i -e 's/>${region}/>${sample}/g' ${sample}_${region}.fasta
    """
}

process SAMTOOLS_STATS {

    label "SAMTOOLS"

    input:
    tuple path(cram), path(crai)
    path(ref_fasta)
    path(ref_fai)

    output:
    path("${cram.simpleName}.cramcov"), emit: cov
    path("${cram.simpleName}.cramstat"), emit: stat

    script:
    """
    samtools coverage --reference ${ref_fasta} ${cram} > ${cram.simpleName}.cramcov
    samtools stats ${cram} > ${cram.simpleName}.cramstat
    """
}
