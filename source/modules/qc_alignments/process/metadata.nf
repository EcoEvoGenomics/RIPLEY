process METADATA_PREPEND_KEY_COLUMN {

    input:
    tuple val(key), path(table)

    output:
    path(table)

    script:
    """
    awk -v key="${key}" 'BEGIN { OFS = "\\t" } { print key, \$0 }' ${table} > tmp && mv tmp ${table}
    """
}
