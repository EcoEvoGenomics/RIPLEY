include { PAIR_CHANNEL_TO_SELF } from "../source/nextflow/workflows/PAIR_CHANNEL_TO_SELF.nf"
include { WRITE_POPULATION_CENSUS; JOIN_GROUPED_CSVS } from "../source/nextflow/modules/system.nf"
include { VCFTOOLS_EXCLUDE_BED } from "../source/nextflow/modules/vcftools.nf"
include { BCFTOOLS_PICK_SAMPLES } from "../source/nextflow/modules/bcftools.nf"
include { REHH_LOAD_VCF; REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY } from "../source/nextflow/modules/rehh.nf"
include { REHH_CALCULATE_IHS; REHH_CALCULATE_XPEHH } from "../source/nextflow/modules/rehh.nf"
include { PLOT_REHH_XPEHH } from "../source/nextflow/modules/plotting.nf"

nextflow.preview.output = true

workflow {
    main:
    pops = Channel.from(params.hs_populations)
    vcfs = Channel.fromPath("${params.hs_vcfdir}/**.vcf.gz")
    vcfs_preprocessed = VCFTOOLS_EXCLUDE_BED(vcfs, params.hs_exclude_bedfile)
    
    pop_vcfs = WRITE_POPULATION_CENSUS(pops, params.metadata) \
    | combine(vcfs_preprocessed) \
    | BCFTOOLS_PICK_SAMPLES

    pop_scans = REHH_LOAD_VCF(pop_vcfs, params.hs_is_polarised) \
    | REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY \
    | flatten \
    | map { scan -> tuple(scan.name.tokenize("_").get(0), scan) } \
    | groupTuple(by: 0) \
    | JOIN_GROUPED_CSVS
    pairwise_pop_scans = PAIR_CHANNEL_TO_SELF(pop_scans)

    REHH_CALCULATE_IHS(
        pop_scans,
        params.hs_is_polarised,
        params.hs_ihs_freqbin,
        params.hs_cand_pval,
        params.hs_cand_window,
        params.hs_cand_overlap,
        params.hs_cand_min_n_mrk,
        params.hs_cand_min_n_extr_mrk,
        params.hs_cand_min_perc_extr_mrk
    )

    REHH_CALCULATE_XPEHH(
        pairwise_pop_scans,
        params.hs_cand_pval,
        params.hs_cand_window,
        params.hs_cand_overlap,
        params.hs_cand_min_n_mrk,
        params.hs_cand_min_n_extr_mrk,
        params.hs_cand_min_perc_extr_mrk
    )

    xpehh_resultfile_list = REHH_CALCULATE_XPEHH.out.csv \
    | map { xpehh -> "${xpehh.toString()}" } \
    | collectFile(name: "xpehh.list", sort: true, newLine: true)

    xpehh_candfile_list = REHH_CALCULATE_XPEHH.out.candidates \
    | map { cand -> "${cand.toString()}" } \
    | collectFile(name: "cand.list", sort: true, newLine: true)

    PLOT_REHH_XPEHH(
        file("${launchDir}/source/R/plot_xpehh.R"),
        xpehh_resultfile_list,
        xpehh_candfile_list,
        params.ref_gff,
        params.ref_chroms_renamed,
        params.hs_cand_pval
    )

    publish:
    haplohh = REHH_LOAD_VCF.out.rds
    pop_scans = JOIN_GROUPED_CSVS.out.joined
    ihs_csv = REHH_CALCULATE_IHS.out.csv
    ihs_rds = REHH_CALCULATE_IHS.out.rds
    ihs_candidate_regions = REHH_CALCULATE_IHS.out.candidates
    xpehh = REHH_CALCULATE_XPEHH.out.csv
    xpehh_candidate_regions = REHH_CALCULATE_XPEHH.out.candidates
    xpehh_parsed_main = PLOT_REHH_XPEHH.output.mainplot
    xpehh_parsed_candplots = PLOT_REHH_XPEHH.output.candplots
    xpehh_parsed_candgenes = PLOT_REHH_XPEHH.output.candgenes
    xpehh_parsed_candregions = PLOT_REHH_XPEHH.output.candregions
}

output {
    haplohh { path "haplotype_selection/chrom_haplohh"}
    pop_scans { path "haplotype_selection/gw_ihh" }
    ihs_csv { path "haplotype_selection/gw_ihs" }
    ihs_rds { path "haplotype_selection/gw_ihs" }
    ihs_candidate_regions { path "haplotype_selection/gw_ihs" }
    xpehh { path "haplotype_selection/gw_xpehh" }
    xpehh_candidate_regions { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_main { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_candplots { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_candgenes { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_candregions { path "haplotype_selection/gw_xpehh" }
}
