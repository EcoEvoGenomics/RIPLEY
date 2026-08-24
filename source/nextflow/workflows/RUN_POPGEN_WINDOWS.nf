include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../modules/genomicsgeneral.nf"
include { BCFTOOLS_LIST_SAMPLES } from "../modules/bcftools.nf"

workflow RUN_POPGEN_WINDOWS {

    take:
    vcfs
    metadata
    focal_populations
    window_size
    step_size
    min_sites

    main:
    repo = GET_GENOMICS_GENERAL()
    geno = GENOMICS_GENERAL_VCF_TO_GENO(repo, vcfs)
    samples = BCFTOOLS_LIST_SAMPLES(vcfs)

    geno_with_sample_list = geno
        .mix(samples)
        .map { file ->
            def key = file.simpleName
            tuple(key, file)
        }
        .groupBy()

    GENOMICS_GENERAL_POPGEN_WINDOWS(
        repo,
        geno_with_sample_list,
        focal_populations.collectFile("pops.txt", newLine: true),
        metadata,
        window_size,
        step_size,
        min_sites,
    )

    emit:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
}
