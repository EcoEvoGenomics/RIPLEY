include { BCFTOOLS_LIST_SAMPLES } from "../../../nf/process/bcftools.nf"
include { METADATA_LIST_POPULATION_MEMBERS } from "../../../nf/process/metadata.nf"
include { alphanumericIssue } from "../../library/filekeys.nf"

workflow PARSE_METADATA {

    take:
    sample_metadata_path
    population_metadata_path
    species_metadata_path
    ploidy_sexes
    focal_population_input
    check_vcf
    check_cram

    main:
    sample_metadata = Channel.fromPath(sample_metadata_path, checkIfExists: true)
    population_metadata = Channel.fromPath(population_metadata_path, checkIfExists: true)
    species_metadata = Channel.fromPath(species_metadata_path, checkIfExists: true)

    samples_in_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[0] }

    populations_in_metadata = population_metadata
        .splitCsv()
        .map { i -> i[0] }

    species_in_metadata = species_metadata
        .splitCsv()
        .map { i -> i[0] }

    unique_samples_in_sample_metadata = samples_in_metadata
        .collect(sort: true)
        .flatten()
        .distinct()
    
    unique_populations_in_sample_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[2] }
        .collect(sort: true)
        .flatten()
        .distinct()

    unique_species_in_sample_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[1] }
        .collect(sort: true)
        .flatten()
        .distinct()

    unique_sexes_in_sample_metadata = sample_metadata
        .splitCsv()
        .map { i -> i[3] }
        .collect(sort: true)
        .flatten()
        .distinct()

    sample_metadata
        .splitCsv()
        .subscribe { entry ->
            if (entry.size() != 4) {
                error("Metadata file ${sample_metadata_path} must have exactly four columns: id, species, population, and sex.")
            }
            def sample_issue = alphanumericIssue(entry[0], "sample", "Metadata file ${sample_metadata_path}")
            if (sample_issue) { error(sample_issue) }
            def species_issue = alphanumericIssue(entry[1], "species", "Metadata file ${sample_metadata_path}")
            if (species_issue) { error(species_issue) }
            def population_issue = alphanumericIssue(entry[2], "population", "Metadata file ${sample_metadata_path}")
            if (population_issue) { error(population_issue) }
            def sex_issue = alphanumericIssue(entry[3], "sex", "Metadata file ${sample_metadata_path}")
            if (sex_issue) { error(sex_issue) }
        }

    population_metadata
        .splitCsv()
        .subscribe { entry ->
            if (entry.size() != 2) {
                error("Metadata file ${population_metadata_path} must have exactly two columns: population and colour.")
            }
            def issue = alphanumericIssue(entry[0], "population", "Metadata file ${population_metadata_path}")
            if (issue) { error(issue) }
            if (!(entry[1] ==~ /^#[0-9A-Fa-f]{6}$/)) {
                error("Metadata file ${population_metadata_path} gives non-hexadecimal colour '${entry[1]}' for population '${entry[0]}'.")
            }
        }

    species_metadata
        .splitCsv()
        .subscribe { entry ->
            if (entry.size() != 2) {
                error("Metadata file ${species_metadata_path} must have exactly two columns: species and colour.")
            }
            def issue = alphanumericIssue(entry[0], "species", "Metadata file ${species_metadata_path}")
            if (issue) { error(issue) }
            if (!(entry[1] ==~ /^#[0-9A-Fa-f]{6}$/)) {
                error("Metadata file ${species_metadata_path} gives non-hexadecimal colour '${entry[1]}' for species '${entry[0]}'.")
            }
        }

    samples_in_metadata
        .collect(sort: true)
        .map { keys -> keys.countBy { key -> key }.findAll { _key, n -> n > 1 }.keySet() as List }
        .subscribe { duplicates ->
            if (duplicates) {
                error("Metadata file ${sample_metadata_path} has duplicate sample entry or entries: ${duplicates.join(', ')}.")
            } 
        }

    populations_in_metadata
        .collect(sort: true)
        .map { keys -> keys.countBy { key -> key }.findAll { _key, n -> n > 1 }.keySet() as List }
        .subscribe { duplicates ->
            if (duplicates) {
                error("Metadata file ${population_metadata_path} has duplicate population entry or entries: ${duplicates.join(', ')}.")
            }
        }

    species_in_metadata
        .collect(sort: true)
        .map { keys -> keys.countBy { key -> key }.findAll { _key, n -> n > 1 }.keySet() as List }
        .subscribe { duplicates ->
            if (duplicates) {
                error("Metadata file ${species_metadata_path} has duplicate species entry or entries: ${duplicates.join(', ')}.")
            }
        }


    // Populations in sample_metadata must be specified in population_metadata
    unique_populations_in_sample_metadata
        .combine(populations_in_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .subscribe { n_uncoloured ->
            if (n_uncoloured > 0) {
                error("Metadata file ${population_metadata_path} lacks entry for ${n_uncoloured} populations.")
            }
        }

    // Species in sample_metadata must be specified in species_metadata
    unique_species_in_sample_metadata
        .combine(species_in_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .subscribe { n_uncoloured ->
            if (n_uncoloured > 0) {
                error("Metadata file ${species_metadata_path} lacks entry for ${n_uncoloured} species.")
            }
        }

    // Silently skipped without a ploidy file because ploidy_sexes then never emits from PARSE_REFERENCE_GENOME
    unique_sexes_in_sample_metadata
        .combine(ploidy_sexes)
        .filter { i -> i[0] !in i[1] }
        .map { i -> i[0] }
        .collect(sort: true)
        .subscribe { absent ->
            if (absent) {
                error("Metadata file ${sample_metadata_path} gives sex value(s) absent from the provided ploidy file: ${absent.unique().join(', ')}.")
            }
        }
    
    // Samples in VCFs must be specified by sample_metadata
    if (check_vcf != null) {
        unique_samples_in_vcf = BCFTOOLS_LIST_SAMPLES(check_vcf)
            .map { sample_list -> sample_list.readLines() }
            .collect(sort: true)
            .flatten()
            .distinct()
        unique_samples_in_vcf
            .combine(unique_samples_in_sample_metadata.toList().toList())
            .filter { i -> i[0] !in i[1] }
            .count()
            .subscribe { n_lacking_metadata ->
                if (n_lacking_metadata > 0) {
                    error("Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
                }
            }
    }

    // Filenames among CRAMs must be specified by sample_metadata
    if (check_cram != null) {
        check_cram.map{ cram -> cram.simpleName }
            .collect()
            .flatten()
            .distinct()
            .combine(unique_samples_in_sample_metadata.toList().toList())
            .filter { i -> i[0] !in i[1] }
            .count()
            .subscribe { n_lacking_metadata ->
                if (n_lacking_metadata > 0) {
                    error("Metadata file ${sample_metadata_path} lacks entry for ${n_lacking_metadata} samples.")
                }
            }
    }

    // All populations are focal unless specified
    focal_populations = (focal_population_input != null)
        ? Channel.from(focal_population_input)
        : unique_populations_in_sample_metadata
    
    // All focal populations must have members in sample_metadata
    focal_populations
        .combine(unique_populations_in_sample_metadata.toList().toList())
        .filter { i -> i[0] !in i[1] }
        .count()
        .subscribe { n_undefined_focal_populations -> 
            if (n_undefined_focal_populations > 0) {
                error("Metadata file ${sample_metadata_path} does not include members for all specified focal populations.")  
            }
        }

    // Census lists for population specification in downstream software
    sample_census = samples_in_metadata.collectFile( name: "samples.list", newLine: true, sort: true )
    focal_populations_censuses = METADATA_LIST_POPULATION_MEMBERS(focal_populations, sample_metadata)

    // Populations maps [id, pop] for file matching in workflows
    sample_population_map = sample_metadata
        .splitCsv()
        .map { entry -> tuple(entry[0], entry[2]) }
    focal_population_map = sample_population_map
        .combine(focal_populations.toList().toList())
        .filter { _sample, population, focal -> population in focal }
        .map { sample, population, _focal -> tuple(sample, population) }

    emit:
    sample_metadata = sample_metadata
    population_metadata = population_metadata
    species_metadata = species_metadata
    sample_census = sample_census
    sample_population_map = sample_population_map
    focal_populations = focal_populations
    focal_populations_censuses = focal_populations_censuses
    focal_population_map = focal_population_map

}
