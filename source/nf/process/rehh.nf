process REHH_LOAD_VCF {

    // NB! Always assumes VCFs are phased but unpolarised.

    label "REHH"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.haplohh.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")

    hh <- rehh::data2haplohh(
        hap_file = "${vcf.toString()}",
        polarize_vcf = FALSE
    )

    saveRDS(hh, file = "${vcf.simpleName}.haplohh.rds")
    """
}

process REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY {

    label "REHH"

    input:
    path(haplohh)

    output:
    path("${haplohh.simpleName}.hh.csv"), emit: csv

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")

    scan <- rehh::scan_hh(
        haplohh = readRDS("${haplohh.toString()}"),
        threads = ${task.cpus}
    )

    write.csv(scan, row.names = FALSE, file = "${haplohh.simpleName}.hh.csv")
    """
}

process REHH_CALCULATE_IHS {

    // NB! Parameter freqbin = 0 assumes original input unpolarised (see process REHH_LOAD_VCF)

    label "REHH"

    input:
    path(csv)
    val(window_size)
    val(step_size)
    val(min_sites)

    output:
    path("${csv.simpleName}.ihs.csv")

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")
    
    ihs <- rehh::ihh2ihs(
        scan = read.csv("${csv.toString()}"),
        freqbin = 0
    )

    windows <- rehh::calc_candidate_regions(
        scan = ihs,
        window_size = ${window_size.toString()},
        overlap = ${step_size.toString()},
        min_n_mrk = ${min_sites.toString()},
        join_neighbors = FALSE,
        threshold = 0
    )

    write.csv(windows, row.names = FALSE, file = "${csv.simpleName}.ihs.csv")
    """
}

process REHH_CALCULATE_XPEHH {

    label "REHH"

    input:
    tuple val(key), val(pop_a), val(pop_b), path(csv_a), path(csv_b)
    val(window_size)
    val(step_size)
    val(min_sites)

    output:
    path("${key}_${pop_a}_${pop_b}.xpehh.csv")

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")
    
    xpehh <- rehh::ies2xpehh(
        scan_pop1 = read.csv("${csv_a.toString()}"),
        scan_pop2 = read.csv("${csv_b.toString()}"),
        popname1 = "${pop_a}",
        popname2 = "${pop_b}",
        include_freq = TRUE
    )

    windows <- rehh::calc_candidate_regions(
        scan = xpehh,
        window_size = ${window_size.toString()},
        overlap = ${step_size.toString()},
        min_n_mrk = ${min_sites.toString()},
        join_neighbors = FALSE,
        threshold = 0
    )

    write.csv(windows, row.names = FALSE, file = "${key}_${pop_a}_${pop_b}.xpehh.csv")
    """
}
