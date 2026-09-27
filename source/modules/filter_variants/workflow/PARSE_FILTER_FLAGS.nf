include { alphanumericIssue; keyFor } from "../../../common/library/filekeys.nf"

workflow PARSE_FILTER_FLAGS {

    take:
    filter_flags_path

    main:
    def permitted_extensions = ["txt"]
    def flags_name = file(filter_flags_path).name

    def flags_key = keyFor(flags_name, permitted_extensions)
    if (flags_key == null) { error("Filter file ${flags_name} must end in .txt.") }
    def issue = alphanumericIssue(flags_key, "filter key", "Filter file ${flags_name}")
    if (issue) { error(issue) }

    filter_flags = Channel.value(file(filter_flags_path, checkIfExists: true))

    emit:
    flags = filter_flags

}
