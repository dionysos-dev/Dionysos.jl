# Inclusion mode for discretizing a continuous set into cells: keep cells fully
# inside (`INNER`), cells that intersect (`OUTER`), or cells whose center lies in
# the set (`CENTER`). Lives in Utils so both `Problem` (AP semantics) and
# `Mapping` (discretization) can reference it despite their load order.

@enum INCL_MODE INNER OUTER CENTER

"Inclusion mode to use for the hole `B` when discretizing `A \\ B` (soundness inverts it)."
@inline function invert_incl_mode(mode::INCL_MODE)
    mode === INNER && return OUTER
    mode === OUTER && return INNER
    mode === CENTER && return OUTER
    return error("Invalid inclusion mode $mode")
end
