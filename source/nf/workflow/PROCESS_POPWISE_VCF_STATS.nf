include { PLOT_POPWISE_VCF_STATS } from "../process/plotting.nf"

workflow PROCESS_POPWISE_VCF_STATS {

    take:
    vcf_stats

    main:
    with_pop = vcf_stats
        .flatten()
        .map { it ->
            def pop = it.simpleName.tokenize("_")[1]
            tuple(pop, it)
        } \
        | PREPEND_POP_COLUMN
    
    across_pop = with_pop
        .flatten()
        .map { it ->
            def ext = it.extension
            def key = it.simpleName.tokenize("_")[0]
            tuple(key, ext, it)
        }
        .collectFile( { it -> ["${it[0]}_popwise.${it[1]}", it[2]] },
            skip: 1,
            keepHeader: true
        )
        .collect()

    plot = PLOT_POPWISE_VCF_STATS(across_pop)

    emit:
    data = across_pop
    plot = plot

}

process PREPEND_POP_COLUMN {

    label "SYSTEM"

    input:
    tuple val(pop), path(table)

    output:
    path(table)

    script:
    """
    sed -e '1s/^/POP\\t/' -e '2,\$s/^/${pop}\\t/' ${table} > tmp && mv tmp ${table} 
    """
}
