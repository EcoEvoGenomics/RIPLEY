include { BEDTOOLS_MAKEWINDOWS } from "../../../common/process/bedtools.nf"
include { RECAST_PIXY_FOR_DCMS; PIXY_STATS } from "../process/pixy.nf"
include { gatedBy } from "../../../common/library/gates.nf"

workflow RUN_PIXY_SCANS {

    take:
    vcfs_indexed
    focal_population_map
    genome_index
    stats
    window_size
    step_size

    main:
    def dcms_contracts = file("${moduleDir}/../library/recast_for_dcms.R", checkIfExists: true)

    def permitted_stats = ["pi", "dxy", "fst", "watterson_theta", "tajima_d"]
    def requested_stats = (stats instanceof List) ? stats : [stats]
    def unknown_stats = requested_stats - permitted_stats
    if (unknown_stats) {
        error("Requested pixy statistic(s) ${unknown_stats.join(', ')} are not among the permitted: ${permitted_stats.join(', ')}.")
    }

    popmap = focal_population_map.collectFile(name: "pixy.popmap.tsv", sort: true) { sample, population ->
        "${sample}\t${population}\n"
    }

    // Population verification guard sits on the dataflow path so pixy cannot start before it lands
    def pairwise_stats = requested_stats.intersect(["dxy", "fst"])
    populations_verified = focal_population_map
        .map { _sample, population -> population }
        .unique()
        .count()
        .map { n_populations ->
            if (pairwise_stats && n_populations < 2) {
                error("Requested pixy statistic(s) ${pairwise_stats.join(', ')} compare populations pairwise, but only ${n_populations} focal population is defined.")
            }
            return true
        }

    windows = BEDTOOLS_MAKEWINDOWS(genome_index, window_size, step_size).bed_base_zero
    pixy_input = vcfs_indexed.combine(windows)
    pixy = PIXY_STATS(pixy_input, gatedBy(popmap, populations_verified).first(), requested_stats)

    // The stat is not carried by the output channel, so it is recovered from the filename
    dcms_formatted = RECAST_PIXY_FOR_DCMS(pixy.map { tsv -> tuple(tsv, tsv.name.tokenize(".")[1]) }, dcms_contracts)

    emit:
    pixy = pixy
    dcms_formatted = dcms_formatted

}
