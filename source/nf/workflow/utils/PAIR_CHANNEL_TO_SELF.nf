workflow PAIR_CHANNEL_TO_SELF {
    
    take:
    input_channel

    main:
    // Pairwise channel self-comparison without item self-comparison by David Mas-Ponte
    // https://github.com/nextflow-io/nextflow/discussions/2109

    pairwise_channel = input_channel
        .combine(input_channel)
        .filter { combined -> combined[0] != combined[1] }
        .map { combined -> combined.sort() }
        .unique()

    emit:
    pairwise_channel

}
