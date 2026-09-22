// Shared helpers for validating the semantic keys RIPLEY encodes into filenames.

// Returns a message rather than raising, because error() from inside a module
// function surfaces only as "Unexpected error [InvocationTargetException]".
def alphanumericIssue(value, kind, source) {
    def text = value?.toString()
    if (text ==~ /^[A-Za-z0-9]+$/) { return null }
    return "${source} gives ${kind} '${text}', which must be strictly alphanumeric (A-Z, a-z, 0-9)."
}

// Returns the name minus its permitted extension, or null if no extension matches.
// Unlike simpleName, which truncates at the first dot, this keeps any interior dot
// in the key so e.g. SAMPLE.A.cram fails validation instead of silently becoming SAMPLE.
def keyFor(name, extensions) {
    def text = name?.toString()
    def matched = extensions.find { ext -> text?.endsWith(".${ext}") }
    if (matched == null) { return null }
    return text[0..<(text.length() - matched.length() - 1)]
}

// Returns the number of tokens in a name when split by a given separator.
def tokenCount(name, separator) {
    def text = name?.toString()
    def tokens = text.tokenize(separator)
    return(tokens.size())
}
