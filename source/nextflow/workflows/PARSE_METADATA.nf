include { BCFTOOLS_LIST_SAMPLES } from "../modules/bcftools.nf"

workflow PARSE_METADATA {

    take:
    sample_metadata_path
    vcfs

    main:
    sample_metadata = Channel.fromPath(sample_metadata_path, checkIfExists: true)
    samples_in_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[0] }
        .toList().toList()

    samples_in_vcfs = BCFTOOLS_LIST_SAMPLES(vcfs)
        .map { sample_list -> sample_list.readLines() }
        .collect(sort:true)
        .flatten()
        .distinct()
        .view()

    samples_in_vcfs
        .combine(samples_in_metadata)
        .filter { i -> i[0] !in i[1] }
        .count()
        .map { n_lacking_metadata ->
            if (n_lacking_metadata > 0) {
                exit(1, "Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
            }
        }

    emit:
    for_samples = sample_metadata

}
