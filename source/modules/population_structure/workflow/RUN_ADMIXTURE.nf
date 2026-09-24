include { PLINK_TO_VCF; PLINK_EXTRACT_SITES } from "../process/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS; CALCULATE_AIM_HIHET } from "../process/admixture.nf"
include { BCFTOOLS_VCF_TO_GENOTABLE } from "../process/bcftools.nf"
include { PLOT_ADMIXTURE; PLOT_HIHET } from "../process/plot.nf"

workflow RUN_ADMIXTURE {

    take:
    plinkfiles
    sample_metadata
    population_metadata
    species_metadata
    kmin
    kmax
    aim_variance_threshold

    main:
    plinkfiles.count().subscribe { n -> n == 1 ?: error("Can only run ADMIXTURE on one VCF.") }
    Channel.of(aim_variance_threshold).subscribe { t ->
        (t as Double) > 0 && (t as Double) < 0.5 ?: error("aim_variance_threshold must be greater than 0 and less than 0.5.")
    }
    
    k_values = Channel.of(kmin..kmax)
    admixture = ADMIXTURE(plinkfiles, k_values)
    admixture_clusts = admixture.clust.collectFile( name: "admixture.clusts" )
    admixture_errors = admixture.error
        .map { cv ->
            def k = cv[0]
            def cv_error = cv[1]
            "${k}\t${cv_error}\n"
        }
        .collectFile( name: "admixture.errors", sort: { line -> line.tokenize("\t")[0] as Integer } )

    numeric_errors = admixture.error
        .filter { cv ->
            def valid_result = (cv[1] as String) ==~ /^-?\d+(\.\d+)?([eE][-+]?\d+)?$/
            valid_result
        }
    best_k = numeric_errors
        .reduce { i, j -> (j[1] as Double) < (i[1] as Double) ? j : i }
        .map { cv -> def k = cv[0]; k }

    admixture_plot = PLOT_ADMIXTURE(
        file("${moduleDir}/../library/plot_admixture.R", checkIfExists: true),
        admixture_clusts, best_k, sample_metadata, population_metadata, species_metadata
    )

    aim_snps = ADMIXTURE_AIMS(admixture.alleles, aim_variance_threshold)
    aim_vcfs = PLINK_EXTRACT_SITES(
        plinkfiles,
        aim_snps.snpids.flatten().filter { snplist -> snplist.readLines().size > 0 }
    ) | PLINK_TO_VCF
    aim_gts = BCFTOOLS_VCF_TO_GENOTABLE(aim_vcfs)

    aim_hihet = aim_snps.alleles
        .flatten()
        .mix(aim_gts)
        .map { it ->
            def stem = it.simpleName.tokenize("_")
            def k = stem[-2][1..-1]
            def pops = stem[-1].tokenize("p")
            def key = "${k}_${pops[0]}_${pops[1]}"
            tuple(key, it)
        }
        .groupTuple(by: 0)
        .map { it -> it[1] }
        .filter { it -> it.size() == 2 } | CALCULATE_AIM_HIHET

    hihet_tables = aim_hihet.hihet
        .map { it ->
            def k = it.simpleName.tokenize("_")[-2][1..-1]
            tuple(k, it)
        }
        .collectFile(
            { it -> ["k${it[0]}.hihet", it[1]] },
            keepHeader: true,
            skip: 1,
            sort: true
        )
        .combine(sample_metadata)
        .combine(population_metadata)

    hihet_plots = PLOT_HIHET(file("${moduleDir}/../library/plot_hihet.R", checkIfExists: true), hihet_tables)

    emit:
    data = admixture.data
    plot = admixture_plot.mix(hihet_plots)
    clusts = admixture_clusts
    errors = admixture_errors
    aims = aim_hihet.aims
    hihet = aim_hihet.hihet

}
