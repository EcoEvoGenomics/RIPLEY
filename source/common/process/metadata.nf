process METADATA_LIST_POPULATION_MEMBERS {

    input:
    each(population)
    path(metadata)
    
    output:
    path("${population}.list")

    script:
    """
    cat ${metadata} | awk -F, '{if (\$3 == "${population}") print \$1}' > "${population}.list"
    if [[ \$(wc -l < "${population}.list") -eq 0 ]]; then
        echo "ERROR: The population ${population} has no members."
        exit 1
    fi
    """
}
