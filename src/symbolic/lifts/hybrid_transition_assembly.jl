# ================================================================
# Transition assembly (pure: operates on already-built per-mode models)
# ================================================================

"""
    build_all_transitions(hs, mode_models, input_mapping) -> (transitions, report)

Assemble every hybrid transition: intra-mode (each mode model's own transitions,
which already encode any time advance) followed by inter-mode (guarded switches
between modes). The [`HybridBuildReport`](@ref) records what the switch wiring
had to discard.
"""
function build_all_transitions(hs::HybridSystem, mode_models, input_mapping::GlobalInputMap)
    transition_list = Vector{HybridTransition}()
    sizehint!(transition_list, sum(get_n_transitions(m) for m in mode_models))

    add_intra_mode_transitions!(transition_list, mode_models, input_mapping)
    report = add_inter_mode_transitions!(transition_list, hs, mode_models, input_mapping)

    return transition_list, report
end

"Embed each mode model's own transitions into the global list, relabelling inputs."
function add_intra_mode_transitions!(
    transition_list,
    mode_models,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(mode_models) "No mode models provided"

    for (mode_id, mode_model) in enumerate(mode_models)
        for (target_local, source_local, local_input_id) in enum_transitions(mode_model)
            global_input_id = get_global_input_id(input_mapping, mode_id, local_input_id)
            global_input_id > 0 || continue
            target_state = (target_local, mode_id)::HybridState
            source_state = (source_local, mode_id)::HybridState
            push!(
                transition_list,
                (target_state, source_state, global_input_id)::HybridTransition,
            )
        end
    end
end

"""
    add_inter_mode_transitions!(transition_list, hs, mode_models, input_mapping) -> HybridBuildReport

Add the guarded mode-switch transitions using each transition's guard and reset map.

A transition that ends up contributing **no** switch is an error, not a warning: the model
declared a switch, the abstraction silently has none, and synthesis then answers for a
disconnected system while reporting success. A *partial* loss is legitimate and is recorded in
the returned report instead.
"""
function add_inter_mode_transitions!(
    transition_list,
    hs::HybridSystem,
    mode_models,
    input_mapping::GlobalInputMap,
)
    transitions = collect(HybridSystems.transitions(hs.automaton))
    report = HybridBuildReport()

    for (transition_id, transition) in enumerate(transitions)
        global_input_id = get_switching_global_id(input_mapping, transition_id)
        global_input_id > 0 || error(
            "Transition $transition_id has no switching input in the global alphabet. The " *
            "input map and the automaton disagree about how many transitions there are.",
        )

        source_mode = HybridSystems.source(hs.automaton, transition)
        target_mode = HybridSystems.target(hs.automaton, transition)

        @assert source_mode <= length(mode_models) "Source mode $source_mode out of bounds"
        @assert target_mode <= length(mode_models) "Target mode $target_mode out of bounds"

        reset_map = HybridSystems.resetmap(hs, transition)
        guard = HybridSystems.guard(hs, transition)

        isnothing(guard) && error(
            "Transition $transition_id (mode $source_mode → mode $target_mode) has no guard, " *
            "so nothing says when the switch may be taken.",
        )

        source_model = mode_models[source_mode]
        target_model = mode_models[target_mode]

        # Source states intersecting the guard (a box over the mode's `[x; t]` coord).
        source_locals = get_states_from_set(source_model, guard, MP.INNER)

        isempty(source_locals) && error(
            "Transition $transition_id (mode $source_mode → mode $target_mode) has a guard no " *
            "cell lies inside, so the switch can never be taken and the modes would be " *
            "disconnected. `INNER` inclusion needs a cell to sit *entirely* within the guard, " *
            "so a guard thinner than one cell — an equality, say — discretizes to nothing. " *
            "Widen the guard into a band, or refine mode $source_mode's `state_grid`.",
        )

        n_dropped = 0
        max_offset = 0.0

        for source_local in source_locals
            # Concretize, apply the reset map, and abstract in the target mode.
            source_coord = get_concrete_state(source_model, source_local)
            reset_result = MS.apply(reset_map, source_coord)
            target_local = abstract_switch_target(target_model, reset_result)

            if target_local <= 0
                n_dropped += 1
                continue
            end

            offset = _reset_quantization_offset(target_model, target_local, reset_result)
            offset === nothing || (max_offset = max(max_offset, offset))

            target_state = (target_local, target_mode)::HybridState
            source_state = (source_local, source_mode)::HybridState
            push!(
                transition_list,
                (target_state, source_state, global_input_id)::HybridTransition,
            )
        end

        if _warn_reset_not_lattice_exact(transition_id, max_offset)
            report.inexact_resets[transition_id] = max_offset
        end

        n_dropped == length(source_locals) && error(
            "Transition $transition_id (mode $source_mode → mode $target_mode): the reset " *
            "image of every one of its $(n_dropped) guard cells falls outside mode " *
            "$target_mode's domain, so the switch can never be taken and the modes would be " *
            "disconnected. Check that the reset lands inside the target mode's bounds.",
        )

        if n_dropped > 0
            report.dropped_resets[transition_id] = (n_dropped, length(source_locals))
            @warn(
                "Transition $transition_id: $n_dropped of $(length(source_locals)) guard " *
                "cells reset outside the target mode's domain; those switches are dropped.",
                maxlog = 1,
            )
        end
    end

    return report
end

# The reset is applied to the source cell *centre* and the image is quantized to
# the nearest target cell, so the abstract switch is exact only when the reset is
# a lattice automorphism — an identity, or a permutation of axes sharing a step,
# or an integer lattice shift. Off the lattice the image is silently snapped by up
# to half a cell, which yields a plausible-looking but unsound abstraction. The
# offset between the image and the centre it snapped to measures exactly that.
function _reset_quantization_offset(m::SymbolicModel, target_local::Int, reset_result)
    target_coord = get_concrete_state(m, target_local)
    scale = max(1.0, maximum(abs, reset_result))
    return maximum(abs, target_coord .- reset_result) / scale
end

# A clock-lifted target rounds the time coordinate *up* by design, so the offset
# is not a soundness signal there.
_reset_quantization_offset(::ClockLiftedSymbolicModel, ::Int, _) = nothing

# Warn if the reset was snapped off-lattice, returning whether it was.
function _warn_reset_not_lattice_exact(transition_id, max_offset; rtol = 1e-9)
    max_offset <= rtol && return false
    @warn(
        "Transition $transition_id: the reset map is not lattice-exact (relative " *
        "offset $(round(max_offset; sigdigits = 3)) between a reset image and the cell " *
        "centre it was quantized to). The switch is applied to cell centres and " *
        "snapped to the nearest cell, so the abstraction may be unsound. Align the " *
        "reset with the grid (identity, permutation of axes sharing a step, or an " *
        "integer lattice shift), or refine the grid.",
        maxlog = 1,
    )
    return true
end

"""
    build_symbolic_automaton(transition_list, mode_models, input_mapping) -> (flat, automaton)

Flatten the hybrid transitions into an [`IndexedAutomatonList`](@ref), returning
`(flat, automaton)` where `flat` is the [`FlatIndex`](@ref) between the integer
numbering and the `(local_id, mode_id)` pairs.
"""
function build_symbolic_automaton(
    transition_list,
    mode_models,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(mode_models) "Mode models cannot be empty"
    # An assertion here reports only that a vector is empty. The cause is always upstream — every
    # mode discretized to nothing, or every guard and every mode transition was dropped — and a
    # model that gets this far deserves to be told which.
    isempty(transition_list) && error(
        "The abstraction has no transitions at all: no mode produced a successor and no switch " *
        "survived its guard. Check that each mode's `state_grid` and `time_step` give the " *
        "dynamics room to move a cell, and that the guards intersect their source mode.",
    )

    estimated_states = sum(get_n_state(m) for m in mode_models)

    states = Set{HybridState}()
    sizehint!(states, estimated_states)

    for (target, source, _) in transition_list
        push!(states, target)
        push!(states, source)
    end

    flat = FlatIndex(collect(states))

    # The automaton is sized by the *alphabet*, not by the inputs that happen to
    # appear: global ids index the whole alphabet, so a mode whose input is never
    # used (a `state_input_filter` pruning it, say) leaves a gap and the largest
    # used id then exceeds a count-based `nsymbols` — the discrete solvers size
    # dense `nstates × nsymbols` tables and index them by symbol.
    symbolic_automaton = IndexedAutomatonList(n_flat(flat), input_mapping.total_inputs)

    @inbounds for (target, source, abstract_input) in transition_list
        target_int = flat_id(flat, target)
        source_int = flat_id(flat, source)
        add_transition!(symbolic_automaton, source_int, target_int, abstract_input)
    end

    return flat, symbolic_automaton
end
