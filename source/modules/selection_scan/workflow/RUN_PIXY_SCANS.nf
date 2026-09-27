include { PIXY_STATS } from "../process/pixy.nf"
include { gatedBy } from "../../../common/library/gates.nf"

workflow RUN_PIXY_SCANS {

    take:
    vcfs_indexed
    focal_population_map
    stats
    window_size

    main:
    def permitted_stats = ["pi", "dxy", "fst", "watterson_theta", "tajima_d"]
    def requested_stats = (stats instanceof List) ? stats : [stats]
    def unknown_stats = requested_stats - permitted_stats
    if (unknown_stats) {
        error("Requested pixy statistic(s) ${unknown_stats.join(', ')} are not among the permitted: ${permitted_stats.join(', ')}.")
    }

    popmap = focal_population_map.collectFile(name: "pixy.popmap.tsv", sort: true) { sample, population ->
        "${sample}\t${population}\n"
    }

    // Guard sits on the dataflow path, not in a subscribe, so pixy cannot start before it lands
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

    PIXY_STATS(vcfs_indexed, gatedBy(popmap, populations_verified).first(), requested_stats, window_size)

    emit:
    pixy = PIXY_STATS.out

}