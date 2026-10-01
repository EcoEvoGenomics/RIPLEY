process CALCULATE_DCMS {

    label "RBASE"

    input:
    path(calculate_dcms_rscript)
    path(tables, stageAs: "dcms/*")

    output:
    path("*.dcms.tsv")

    script:
    """
    Rscript ${calculate_dcms_rscript}
    """
}
