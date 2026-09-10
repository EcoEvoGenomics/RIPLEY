include { BCFTOOLS_LIST_SAMPLES } from "../process/bcftools.nf"

workflow PARSE_METADATA {

    // To-do: Add test to ensure all metadata is strictly alphanumeric

    take:
    sample_metadata_path
    focal_population_input
    check_vcf
    check_cram

    main:
    sample_metadata = Channel.fromPath(sample_metadata_path, checkIfExists: true)

    // Samples must have no duplicate entries in metadata
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

    // If user set no focal populations, all are focal
    unique_populations_in_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[2] }
        .collect(sort: true)
        .flatten()
        .distinct()

    focal_populations_set_by_user = (focal_population_input != null)
        ? Channel.of(true)
        : Channel.of(false)
        
    focal_populations = (focal_population_input != null)
        ? Channel.from(focal_population_input)
        : unique_populations_in_metadata
    
    // All focal populations must have members in metadata
    focal_populations
        .combine(unique_populations_in_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .map { n_undefined_focal_populations -> 
            if (n_undefined_focal_populations > 0) {
                exit(1, "Metadata file ${sample_metadata_path} does not include members for all specified focal populations.")  
            }
        }
    
    // All samples in VCF channel must be specified by metadata
    if (check_vcf != null) {
        unique_samples_in_vcf = BCFTOOLS_LIST_SAMPLES(check_vcf)
            .map { sample_list -> sample_list.readLines() }
            .collect(sort: true)
            .flatten()
            .distinct()
        unique_samples_in_vcf
            .combine(unique_samples_in_metadata.toList().toList())
            .filter { i -> i[0] !in i[1] }
            .count()
            .map { n_lacking_metadata ->
                if (n_lacking_metadata > 0) {
                    exit(1, "Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
                }
            }
    }

    // All filenames in CRAM channel must be specified by metadata
    if (check_cram != null) {
        check_cram.map{ cram -> cram.simpleName }
            .collect()
            .flatten()
            .distinct()
            .combine(unique_samples_in_metadata.toList().toList())
            .filter { i -> i[0] !in i[1] }
            .count()
            map { n_lacking_metadata ->
                if (n_lacking_metadata > 0) {
                    exit(1, "Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
                }
            }
    }
    
    emit:
    focal_populations_set_by_user = focal_populations_set_by_user
    focal_populations = focal_populations
    for_samples = sample_metadata

}
