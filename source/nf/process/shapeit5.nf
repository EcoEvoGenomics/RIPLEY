process SHAPEIT5_PHASE_COMMON {

    label "SHAPEIT5"

    input:
    tuple val(name), val(region), path(vcf), path(csi)

    output:
    tuple path("${name}.bcf"), path("${name}.bcf.csi")

    script:
    """
    SHAPEIT5_phase_common \
        --input ${vcf} \
        --region ${region} \
        --output ${name}.bcf \
        --thread ${task.cpus}
    """
}

process SHAPEIT5_LIGATE {

    // Window keys are zero-padded, so lexical order is window order

    label "SHAPEIT5"

    input:
    tuple val(chrom), path(chunks, stageAs: "chunks/*")

    output:
    tuple path("${chrom}.bcf"), path("${chrom}.bcf.csi")

    script:
    """
    ls -1 chunks/*.bcf | sort > chunks.txt

    # SHAPEIT5_ligate requires at least two files; single chunk passes through
    if [ \$(wc -l < chunks.txt) -eq 1 ]
    then
        cp -L \$(cat chunks.txt) ${chrom}.bcf
        cp -L \$(cat chunks.txt).csi ${chrom}.bcf.csi
    else
        SHAPEIT5_ligate \
            --input chunks.txt \
            --output ${chrom}.bcf \
            --index \
            --thread ${task.cpus}
    fi
    """
}
