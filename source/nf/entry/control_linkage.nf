include { PLINK_INIT_PLINKFILES; PLINK_PAIRWISE_LD; PARSE_PLINK_LD_DECAY } from "../process/plink.nf"
include { PLOT_PLINK_LD_DECAY } from "../process/plot.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.cl_vcf)
    
    PLINK_INIT_PLINKFILES(vcf, params.ref_n_chroms)
    PLINK_PAIRWISE_LD(PLINK_INIT_PLINKFILES.out, params.cl_thin, params.cl_window, params.cl_window_kb)
    PARSE_PLINK_LD_DECAY(PLINK_PAIRWISE_LD.out, params.ref_exclude_prefix, params.cl_bin_size)
    PLOT_PLINK_LD_DECAY(
        file("${moduleDir}/../../R/plot_linkage.R", checkIfExists: true),
        PARSE_PLINK_LD_DECAY.out
    )

    publish:
    stats = PLINK_PAIRWISE_LD.out
    decay = PARSE_PLINK_LD_DECAY.out
    decay_plot = PLOT_PLINK_LD_DECAY.out
}

output {
    stats { path "control_linkage" }
    decay { path "control_linkage" }
    decay_plot { path "control_linkage" }
}
