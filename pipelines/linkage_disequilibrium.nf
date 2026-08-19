include { PLINK_INIT_BEDFILES; PLINK_PAIRWISE_LD; PARSE_PLINK_LD_DECAY } from "../source/nextflow/modules/plink.nf"
include { PLOT_PLINK_LD_DECAY } from "../source/nextflow/modules/plotting.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.ld_vcf)
    
    PLINK_INIT_BEDFILES(vcf, params.ref_n_chroms)
    PLINK_PAIRWISE_LD(PLINK_INIT_BEDFILES.out, params.ld_thin, params.ld_window, params.ld_window_kb)
    PARSE_PLINK_LD_DECAY(PLINK_PAIRWISE_LD.out, params.ref_exclude_prefix, params.ld_bin_size)
    PLOT_PLINK_LD_DECAY(PARSE_PLINK_LD_DECAY.out)

    publish:
    ld_stats = PLINK_PAIRWISE_LD.out
    ld_decay = PARSE_PLINK_LD_DECAY.out
    ld_decay_plot = PLOT_PLINK_LD_DECAY.out
}

output {
    ld_stats { path "linkage_disequilibrium" }
    ld_decay { path "linkage_disequilibrium" }
    ld_decay_plot { path "linkage_disequilibrium" }
}
