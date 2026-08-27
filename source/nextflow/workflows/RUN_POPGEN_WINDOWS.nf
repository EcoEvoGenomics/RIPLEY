include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../modules/genomicsgeneral.nf"
include { BCFTOOLS_LIST_SAMPLES } from "../modules/bcftools.nf"

workflow RUN_POPGEN_WINDOWS {

    take:
    vcfs
    metadata
    populations
    window_size
    step_size
    min_sites

    main:
    repo = GET_GENOMICS_GENERAL()
    geno = GENOMICS_GENERAL_VCF_TO_GENO(repo, vcfs)
    samples = BCFTOOLS_LIST_SAMPLES(vcfs)

    geno_with_sample_list = geno
        .combine(samples)
        .map { files ->
            def key_geno = files[0].simpleName
            def key_samples = files[1].simpleName
            key_geno == key_samples
                ? tuple(files[0], files[1])
                : null
        }

    GENOMICS_GENERAL_POPGEN_WINDOWS(
        repo,
        geno_with_sample_list,
        populations.collectFile(name: "pops.txt", newLine: true),
        metadata,
        window_size,
        step_size,
        min_sites,
    )

    emit:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
}
