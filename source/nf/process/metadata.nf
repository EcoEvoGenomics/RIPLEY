process METADATA_PREPEND_KEY_COLUMN {

    label "BASE"

    input:
    tuple val(key), path(table)

    output:
    path(table)

    script:
    """
    awk -v key="${key}" 'BEGIN { OFS = "\\t" } { print key, \$0 }' ${table} > tmp && mv tmp ${table}
    """
}

process METADATA_PREPEND_KEY_COLUMN_WITH_HEADER {

    label "BASE"

    input:
    tuple val(header), val(key), path(table)

    output:
    path(table)

    script:
    """
    awk -v header="${header}" -v key="${key}" 'BEGIN { OFS = "\\t" } NR == 1 { print header, \$0; next } { print key, \$0 }' ${table} > tmp && mv tmp ${table}
    """
}

process METADATA_LIST_POPULATION_MEMBERS {

    label "BASE"

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
