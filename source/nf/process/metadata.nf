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
