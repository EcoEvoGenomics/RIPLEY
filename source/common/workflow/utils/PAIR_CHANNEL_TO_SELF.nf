workflow PAIR_CHANNEL_TO_SELF {

    // TAKE --------------------------------------------------------
    // input_channel    <- [A, B, C]
    // EMIT --------------------------------------------------------
    // pairwise_channel -> [[A, B], [A, C], [B, C]]
    // -------------------------------------------------------------
    //
    // Pairwise channel self-comparison without item self-comparison
    // by David Mas-Ponte. See GitHub:
    // https://github.com/nextflow-io/nextflow/discussions/2109

    take:
    input_channel

    main:
    pairwise_channel = input_channel
        .combine(input_channel)
        .filter { combined -> combined[0] != combined[1] }
        .map { combined -> combined.sort() }
        .unique()

    emit:
    pairwise_channel

}
