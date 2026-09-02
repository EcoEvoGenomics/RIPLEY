include { PARSE_REFERENCE_GENOME } from "../source/nextflow/workflows/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../source/nextflow/workflows/PARSE_VCF.nf"
include { PARSE_METADATA } from "../source/nextflow/workflows/PARSE_METADATA.nf"
include { RUN_KINSHIP_ANALYSIS } from "../source/nextflow/workflows/RUN_KINSHIP_ANALYSIS.nf"
include { RUN_PAIRWISE_FST } from "../source/nextflow/workflows/RUN_PAIRWISE_FST.nf"
include { RUN_LD_PRUNING } from "../source/nextflow/workflows/RUN_LD_PRUNING.nf"
include { RUN_ADMIXTURE } from "../source/nextflow/workflows/RUN_ADMIXTURE.nf"
include { RUN_PCA } from "../source/nextflow/workflows/RUN_PCA.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.ps_vcf, params.ref_exclude_coords, genome.chrom_names, false)
    metadata = PARSE_METADATA(params.metadata, params.focal_populations, input.vcf_condensed)
    pruned = RUN_LD_PRUNING(input.as_plinkfiles, params.ps_prune_window_kb, params.ps_prune_step_snps, params.ps_prune_threshold)

    RUN_KINSHIP_ANALYSIS(input.vcf_condensed, metadata.for_samples)
    RUN_PAIRWISE_FST(input.vcf_condensed, metadata.focal_populations, metadata.for_samples)

    RUN_PCA(pruned.plinkfiles, metadata.for_samples)
    RUN_ADMIXTURE(pruned.plinkfiles, metadata.for_samples, params.ps_admixture_kmin, params.ps_admixture_kmax, params.ps_aim_variance_threshold)

    publish:
    kinship_data = RUN_KINSHIP_ANALYSIS.out.data
    kinship_plot = RUN_KINSHIP_ANALYSIS.out.plot
    fst_data = RUN_PAIRWISE_FST.out.data
    fst_logs = RUN_PAIRWISE_FST.out.logfile
    fst_mean = RUN_PAIRWISE_FST.out.mean
    fst_plot = RUN_PAIRWISE_FST.out.plot
    pca_data = RUN_PCA.out.data
    pca_plot = RUN_PCA.out.plot
    admixture_data = RUN_ADMIXTURE.out.data
    admixture_aims = RUN_ADMIXTURE.out.aims
    admixture_clusts = RUN_ADMIXTURE.out.clusts
    admixture_errors = RUN_ADMIXTURE.out.errors

}

output {

    kinship_data { path "population_structure/kinship" }
    kinship_plot { path "population_structure/kinship" }
    fst_data { path "population_structure/fst/data" }
    fst_logs { path "population_structure/fst/data" }
    fst_mean { path "population_structure/fst" }
    fst_plot { path "population_structure/fst" }
    pca_data { path "population_structure/pca" }
    pca_plot { path "population_structure/pca" }
    admixture_data { path "population_structure/admixture/data" }
    admixture_aims { path "population_structure/admixture/aims" }
    admixture_errors { path "population_structure/admixture" }
    admixture_clusts { path "population_structure/admixture" }

}
