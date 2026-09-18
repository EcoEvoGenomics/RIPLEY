// Shared helpers for validating the semantic keys RIPLEY encodes into filenames.

// Keys are recovered from filenames by splitting on underscores and dots, so a
// key containing either (or whitespace, or a shell metacharacter) is corrupted.
// Returns a message rather than raising, because error() from inside a module
// function surfaces only as "Unexpected error [InvocationTargetException]".
def alphanumericIssue(value, kind, source) {
    def text = value?.toString()
    if (text ==~ /^[A-Za-z0-9]+$/) { return null }
    return "${source} gives ${kind} '${text}', which must be strictly alphanumeric (A-Z, a-z, 0-9)."
}
