"""
    ModeLift{H, G} <: AbstractLift

The mode lift: composes one `AbstractSymbolicModel` per mode into a single
[`HybridSymbolicModel`](@ref), wiring the guarded, reset-coupled switches of the
hybrid automaton `hs` and unifying inputs through `input_mapping`.

Unlike [`ClockLift`](@ref) — which lifts a *single* base model — `ModeLift`
consumes a *vector* of per-mode base models (each plain or clock-lifted), so it is
applied via `lift(l, mode_models)`.
"""
struct ModeLift{H <: HybridSystem, G <: GlobalInputMap} <: AbstractLift
    hs::H
    input_mapping::G
end

"""
    lift(l::ModeLift, mode_models) -> HybridSymbolicModel

Compose the per-mode `mode_models` (each plain or clock-lifted) into a hybrid
symbolic model, assembling the intra-mode and guarded inter-mode transitions.
"""
function lift(l::ModeLift, mode_models::AbstractVector)
    transition_list, report = build_all_transitions(l.hs, mode_models, l.input_mapping)
    flat, automaton =
        build_symbolic_automaton(transition_list, mode_models, l.input_mapping)
    return HybridSymbolicModel(mode_models, flat, automaton, l.input_mapping, report)
end
