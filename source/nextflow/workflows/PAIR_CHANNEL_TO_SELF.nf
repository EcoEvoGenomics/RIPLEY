workflow PAIR_CHANNEL_TO_SELF {
    
    take:
    input_channel

    main:
    // Pairwise channel self-comparison without item self-comparison by David Mas-Ponte
    // https://github.com/nextflow-io/nextflow/discussions/2109

    pairwise_channel = input_channel
        .combine(input_channel)
        .filter { scan -> scan[0] != scan[1] }
        .map { scan -> scan.sort() }
        .unique()

    emit:
    pairwise_channel

}
