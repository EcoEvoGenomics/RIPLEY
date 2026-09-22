// Shared helpers for normalising pipeline inputs.

// Coerce a parameter that may be given as null, a scalar, or a list into a list.
def asList(value) {
    if (value == null) { return [] }
    if (value instanceof List) { return value as List }
    return [value]
}
