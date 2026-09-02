include { BCFTOOLS_LIST_SAMPLES } from "../process/bcftools.nf"

workflow PARSE_METADATA {

    take:
    sample_metadata_path
    focal_population_input
    vcfs

    main:
    sample_metadata = Channel.fromPath(sample_metadata_path, checkIfExists: true)

    // To-do: Add test to ensure all metadata is strictly alphanumeric

    // Raise error if there are duplicate sample entries in metadata
    samples_in_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[0] }

    unique_samples_in_metadata = samples_in_metadata
        .collect(sort: true)
        .flatten()
        .distinct()

    unique_samples_in_metadata.count()
        .combine(samples_in_metadata.count())
        .map { counts ->
            if (counts[0] != counts[1]) {
                exit(1, "Metadata file ${sample_metadata_path} has duplicate sample entry or entries.")
            } 
        }

    // Raise error if not all samples in VCFs are in metadata
    unique_samples_in_vcfs = BCFTOOLS_LIST_SAMPLES(vcfs)
        .map { sample_list -> sample_list.readLines() }
        .collect(sort: true)
        .flatten()
        .distinct()

    unique_samples_in_vcfs
        .combine(unique_samples_in_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .map { n_lacking_metadata ->
            if (n_lacking_metadata > 0) {
                exit(1, "Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
            }
        }

    // Parse focal populations: if none provided, use all in metadata
    unique_populations_in_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[2] }
        .collect(sort: true)
        .flatten()
        .distinct()
        
    focal_populations = (focal_population_input != null)
        ? Channel.from(focal_population_input)
        : unique_populations_in_metadata
    
    // Raise error if not all focal populations have members in metadata
    focal_populations
        .combine(unique_populations_in_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .map { n_undefined_focal_populations -> 
            if (n_undefined_focal_populations > 0) {
                exit(1, "Metadata file ${sample_metadata_path} does not include members for all specified focal populations.")  
            }
        }
    
    emit:
    focal_populations = focal_populations
    for_samples = sample_metadata

}
