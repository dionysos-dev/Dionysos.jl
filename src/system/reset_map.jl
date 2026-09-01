
# ------------------------------
# Hybrid transitions: guard + reset
# ------------------------------

"""
    GuardedResetMap{G, F} <: MathematicalSystems.AbstractMap

The map taken on a hybrid transition: it is enabled on the **guard** set `G` and applies the
**reset** `F` to the state.

`MathematicalSystems.stateset` returns the guard — which is how `HybridSystems` stores the
enabling condition of a transition — and `MathematicalSystems.apply` runs the reset. `reset`
defaults to the identity, the common case where switching mode leaves the state untouched.

The state it is applied to is the *augmented* state of the mode: `x` for a plain mode, and
`[x; t]` for a clock-lifted one. `guard` must live in the same space.

```julia
# switch enabled on T ≤ 19, leaving the temperature unchanged
GuardedResetMap(LazySets.Hyperrectangle(; low = [17.0], high = [19.0]))

# switch enabled on the target box, resetting x and clamping the clock
GuardedResetMap(target, state -> vcat([0.0], max(1.0, state[end])))
```
"""
struct GuardedResetMap{G, F} <: MS.AbstractMap
    guard::G
    reset::F
end

GuardedResetMap(guard) = GuardedResetMap(guard, identity)

MS.stateset(map::GuardedResetMap) = map.guard
MS.apply(map::GuardedResetMap, state::AbstractVector) = map.reset(state)
