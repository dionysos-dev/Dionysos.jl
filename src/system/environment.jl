# Who owns which input signal. The ecosystem already says it — `MathematicalSystems` through the
# `isnoisy` trait and `noiseset`, `HybridSystems` through the switching type of Liberzon §1.1.3 —
# and these accessors turn those vocabularies into the one question the rest of the toolbox asks.
# Nothing downstream mentions `isnoisy` or `AutonomousSwitching`; extending to a new ownership
# vocabulary is one more method here.

"""
    environment_input(system)

The environment's choice set of `system` — the set from which an adversary picks at every step —
or `nothing` when every input is the controller's.

This is the single seam between how a system *declares* ownership and how the toolbox *acts* on
it: an abstraction that receives a system with a non-`nothing` environment folds those choices
into its successor sets, and the unchanged solvers then compute the robust game `∃u ∀w`.

- A `MathematicalSystems` system owns an environment iff it is noisy; the set is `MS.noiseset`,
  so only the `Constrained…` variants qualify (the toolbox needs the set, not just the dimension).
- A `HybridSystems.HybridSystem` owns one iff its switching is autonomous; the set is the mode
  labels `1:n`. Mixed ownership across transitions is refused rather than guessed.
- Anything else declares nothing and defaults to `nothing` — pure synthesis, exactly as today.
"""
environment_input(system) = nothing

function environment_input(system::MS.AbstractSystem)
    MS.isnoisy(system) === true || return nothing
    hasmethod(MS.noiseset, Tuple{typeof(system)}) || error(
        "$(typeof(system).name.name) is noisy but carries no noise set; robust synthesis needs " *
        "the set, so use the Constrained… variant of the type.",
    )
    return MS.noiseset(system)
end

function environment_input(system::HybridSystems.HybridSystem)
    switchings = system.switchings
    if all(sw -> sw isa HybridSystems.AutonomousSwitching, switchings)
        return 1:length(system.resetmaps)
    elseif all(sw -> sw isa HybridSystems.ControlledSwitching, switchings)
        return nothing
    end
    return error(
        "The switching ownership is mixed across transitions; fold or enumerate cannot be " *
        "decided for the system as a whole.",
    )
end

"""
    with_switching(system::HybridSystems.HybridSystem, switching)

The same hybrid system with every switching declared as `switching` — the ownership of the mode
signal made explicit.

`HybridSystems.discreteswitchedsystem` hardcodes `AutonomousSwitching`, under which the mode is
the environment's and a solve is a verification; a switched system whose modes the *controller*
picks must say so:

    f = with_switching(HybridSystems.discreteswitchedsystem(A), HybridSystems.ControlledSwitching())

This is the interim spelling until `discreteswitchedsystem` accepts the switching upstream.
"""
function with_switching(
    system::HybridSystems.HybridSystem,
    switching::HybridSystems.AbstractSwitching,
)
    return HybridSystems.HybridSystem(
        system.automaton,
        system.modes,
        system.resetmaps,
        FillArrays.Fill(switching, length(system.switchings)),
        system.ext,
    )
end
