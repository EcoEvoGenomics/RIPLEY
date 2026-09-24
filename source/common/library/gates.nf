// Shared helpers for gating workflow emits on conditions

// Withholds every emission until the gate emits, so no consumer can run ahead of the check.
// The wrapping list keeps tuple emissions intact, which combine would otherwise flatten.
def gatedBy(channel, gate) {
    channel.map { item -> [item] }.combine(gate).map { pair -> pair[0] }
}
