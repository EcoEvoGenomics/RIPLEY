include { PARSE_REFERENCE_GENOME } from "../source/nextflow/workflows/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF_FILE } from "../source/nextflow/workflows/PARSE_VCF_FILE.nf"
include { RUN_KINSHIP_ANALYSIS } from "../source/nextflow/workflows/RUN_KINSHIP_ANALYSIS.nf"
include { RUN_PAIRWISE_FST } from "../source/nextflow/workflows/RUN_PAIRWISE_FST.nf"
include { RUN_LD_PRUNING } from "../source/nextflow/workflows/RUN_LD_PRUNING.nf"
include { RUN_ADMIXTURE } from "../source/nextflow/workflows/RUN_ADMIXTURE.nf"
include { RUN_PCA } from "../source/nextflow/workflows/RUN_PCA.nf"

nextflow.preview.output = true

workflow {
    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix)
    input = PARSE_VCF_FILE(params.ps_vcf, genome.n_chroms, genome.excluded)
    vcf = input.vcf
    bed = input.bedfiles

    // UNPRUNED ANALYSES
    RUN_KINSHIP_ANALYSIS(vcf)
    RUN_PAIRWISE_FST(vcf, params.ps_fst_populations, params.metadata)

    // LD-PRUNED ANALYSES
    pruned = RUN_LD_PRUNING(bed, params.ps_prune_window_kb, params.ps_prune_step_snps, params.ps_prune_threshold).bedfiles
    RUN_ADMIXTURE(pruned, params.ps_admixture_kmin, params.ps_admixture_kmax, params.ps_aim_variance_threshold)
    RUN_PCA(pruned)

    publish:
    kinship_table = RUN_KINSHIP_ANALYSIS.out.data
    kinship_matrix = RUN_KINSHIP_ANALYSIS.out.plot
    pairwise_mean_fst = RUN_PAIRWISE_FST.out.mean
    admixture = RUN_ADMIXTURE.out.data
    aims = RUN_ADMIXTURE.out.aims
    pca = RUN_PCA.out
}

output {
    kinship_table { path "population_structure/kinship" }
    kinship_matrix { path "population_structure/kinship" }
    pairwise_mean_fst { path "population_structure/fst" }
    admixture { path "population_structure/admixture" }
    aims { path "population_structure/aims" }
    pca { path "population_structure/pca" }
}
