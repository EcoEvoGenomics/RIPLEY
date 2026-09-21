include { PARSE_REFERENCE_GENOME } from "../workflow/parse/PARSE_REFERENCE_GENOME.nf"
include { PARSE_VCF } from "../workflow/parse/PARSE_VCF.nf"
include { PARSE_VCF_TO_PLINK } from "../workflow/parse/PARSE_VCF_TO_PLINK.nf"
include { PARSE_METADATA } from "../workflow/parse/PARSE_METADATA.nf"
include { RUN_KINSHIP_ANALYSIS } from "../workflow/run/RUN_KINSHIP_ANALYSIS.nf"
include { RUN_PAIRWISE_FST } from "../workflow/run/RUN_PAIRWISE_FST.nf"
include { RUN_LD_PRUNING } from "../workflow/run/RUN_LD_PRUNING.nf"
include { RUN_ADMIXTURE } from "../workflow/run/RUN_ADMIXTURE.nf"
include { RUN_PCA } from "../workflow/run/RUN_PCA.nf"

nextflow.preview.output = true

workflow {

    main:
    genome = PARSE_REFERENCE_GENOME(params.ref_genome, params.ref_ploidy, params.ref_exclude_chroms, params.ref_exclude_prefix, params.ref_chrom_labels)
    input = PARSE_VCF(params.ps_vcf, params.ref_exclude_coords, genome.chrom_names, true, false)
    metadata = PARSE_METADATA(params.sample_metadata, params.population_metadata, params.species_metadata, genome.ploidy_sexes, params.focal_populations, input.vcf, null)

    plink = PARSE_VCF_TO_PLINK(input.vcf, genome.total_chroms)
    pruned = RUN_LD_PRUNING(plink.as_plinkfiles, params.ps_prune_window_kb, params.ps_prune_step_snps, params.ps_prune_threshold)

    RUN_KINSHIP_ANALYSIS(plink.vcf_condensed, metadata.sample_metadata, metadata.population_metadata, metadata.species_metadata)
    RUN_PAIRWISE_FST(plink.vcf_condensed, metadata.focal_populations_censuses)
    RUN_PCA(pruned.plinkfiles, metadata.sample_metadata, metadata.population_metadata, metadata.species_metadata)
    RUN_ADMIXTURE(
        pruned.plinkfiles, 
        metadata.sample_metadata, 
        metadata.population_metadata,
        metadata.species_metadata,
        params.ps_admixture_kmin, 
        params.ps_admixture_kmax, 
        params.ps_aim_variance_threshold
    )

    publish:
    kinship_data = RUN_KINSHIP_ANALYSIS.out.data
    kinship_plot = RUN_KINSHIP_ANALYSIS.out.plot
    fst_data = RUN_PAIRWISE_FST.out.data
    fst_logs = RUN_PAIRWISE_FST.out.logs
    fst_mean = RUN_PAIRWISE_FST.out.mean
    fst_plot = RUN_PAIRWISE_FST.out.plot
    pca_data = RUN_PCA.out.data
    pca_logs = RUN_PCA.out.logs
    pca_plot = RUN_PCA.out.plot
    admixture_plot = RUN_ADMIXTURE.out.plot
    admixture_data = RUN_ADMIXTURE.out.data
    admixture_aims = RUN_ADMIXTURE.out.aims
    admixture_hihet = RUN_ADMIXTURE.out.hihet
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
    pca_data { path "population_structure/pca/data" }
    pca_logs { path "population_structure/pca" }
    pca_plot { path "population_structure/pca" }
    admixture_plot { path "population_structure/admixture" }
    admixture_data { path "population_structure/admixture/data" }
    admixture_aims { path "population_structure/admixture/aims" }
    admixture_hihet { path "population_structure/admixture/aims" }
    admixture_errors { path "population_structure/admixture" }
    admixture_clusts { path "population_structure/admixture" }

}
