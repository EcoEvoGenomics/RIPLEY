include { BCFTOOLS_LIST_SAMPLES } from "../../process/bcftools.nf"
include { METADATA_LIST_POPULATION_MEMBERS } from "../../process/metadata.nf"
include { alphanumericIssue } from "../../library/filekeys.nf"
include { gatedBy } from "../../library/gates.nf"

// Returns a message rather than raising, because error() from inside a function
// surfaces only as "Unexpected error [InvocationTargetException]".
def reconciliationIssue(in_input, in_metadata, source, kind) {
    def lacking_metadata = in_input - in_metadata
    if (lacking_metadata) {
        return "Metadata file ${source} lacks entry for ${lacking_metadata.size()} samples."
    }
    def absent = in_metadata - in_input
    if (absent) {
        return "Metadata file ${source} gives sample(s) absent from the provided ${kind}(s): ${absent.join(', ')}."
    }
    return null
}

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
    
    samples_listed_in_metadata = unique_samples_in_sample_metadata
        .collect(sort: true)
        .map { samples -> [samples] }

    // Reconciliation sits on the dataflow path, not in a subscribe, so that gated
    // consumers cannot complete before a mismatch is detected
    if (check_vcf != null) {
        samples_verified = BCFTOOLS_LIST_SAMPLES(check_vcf)
            .map { sample_list -> sample_list.readLines() }
            .collect(sort: true)
            .map { samples -> [samples.flatten().unique().sort()] }
            .combine(samples_listed_in_metadata)
            .map { in_vcf, in_metadata ->
                def issue = reconciliationIssue(in_vcf, in_metadata, sample_metadata_path, "VCF")
                if (issue) { error(issue) }
                return true
            }
            .first()
    } else if (check_cram != null) {
        samples_verified = check_cram
            .map { cram -> cram.simpleName }
            .collect(sort: true)
            .map { samples -> [samples.unique().sort()] }
            .combine(samples_listed_in_metadata)
            .map { in_cram, in_metadata ->
                def issue = reconciliationIssue(in_cram, in_metadata, sample_metadata_path, "CRAM")
                if (issue) { error(issue) }
                return true
            }
            .first()
    } else {
        samples_verified = Channel.value(true)
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
    focal_populations_censuses = METADATA_LIST_POPULATION_MEMBERS(
        gatedBy(focal_populations, samples_verified),
        gatedBy(sample_metadata, samples_verified)
    )

    // Populations maps [id, pop] for file matching in workflows
    sample_population_map = sample_metadata
        .splitCsv()
        .map { entry -> tuple(entry[0], entry[2]) }
    focal_population_map = sample_population_map
        .combine(focal_populations.toList().toList())
        .filter { _sample, population, focal -> population in focal }
        .map { sample, population, _focal -> tuple(sample, population) }

    emit:
    sample_metadata = gatedBy(sample_metadata, samples_verified)
    population_metadata = gatedBy(population_metadata, samples_verified)
    species_metadata = gatedBy(species_metadata, samples_verified)
    sample_census = gatedBy(sample_census, samples_verified)
    sample_population_map = gatedBy(sample_population_map, samples_verified)
    focal_populations = gatedBy(focal_populations, samples_verified)
    focal_populations_censuses = focal_populations_censuses
    focal_population_map = gatedBy(focal_population_map, samples_verified)

}
