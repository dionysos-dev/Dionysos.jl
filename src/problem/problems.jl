"""
    ProblemType

An abstract type that represents a generic control problem.  
All concrete problem types (e.g., [`EmptyProblem`](@ref), [`OptimalControlProblem`](@ref), [`SafetyProblem`](@ref), [`CoSafeLTLProblem`](@ref)) should subtype `ProblemType`.
"""
abstract type ProblemType end

"""
    EmptyProblem{S, X} <: ProblemType

A problem type used to construct an **abstraction** of a dynamical system.

- `S`: The system to abstract (continuous or discrete-time).
- `X`: The region of interest (e.g., a subset of the state space).

This problem encodes no control objective. It is intended for generating symbolic models that can later be reused by other solvers.
"""
mutable struct EmptyProblem{S, X} <: ProblemType
    system::S
    region::X
end

"""
    OptimalControlProblem{S, XI, XT, XC, TC, T} <: ProblemType

Encodes a **reach-avoid optimal control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached.
- `XC`: A state cost function or structure.
- `TC`: A transition cost function or structure.
- `T`: The time horizon (number of allowed time steps).

This problem aims to find a control strategy that reaches the target set from the initial set, minimizing the accumulated cost over time.
"""
mutable struct OptimalControlProblem{S, XI, XT, XC, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::Any
    time::T
end

"""
    SafetyProblem{S, XI, XS, T} <: ProblemType

Encodes a **safety control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XS`: The safe set in which the system must remain.
- `T`: The time horizon (number of allowed time steps).

This problem aims to synthesize a controller that ensures the system remains within the safe set for the entire duration of the time horizon.
"""
mutable struct SafetyProblem{S, XI, XS, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

"""
    CoSafeLTLProblem{S, XI, SPEC, LAB} <: ProblemType

Encodes a **co-safe LTL control problem**.

- `S`: The system to control.
- `XI`: The initial set of states.
- `SPEC`: The co-safe LTL specification object (e.g. a Spot/LTL formula or an automaton/monitor wrapper).
- `LAB`: The labeling payload type used in `labeling` (typically a concrete set type such as a LazySet,
         or an abstract labeling such as `Vector{Int}` / bitset / etc.).

# Fields
- `system::S`:
  The (concrete or abstract) system to control.

- `initial_set::XI`:
  Initial set of states (or initial abstract states).

- `spec::SPEC`:
  The co-safe LTL specification.

- `labeling::Dict{Symbol, LAB}`:
  Unified container mapping each atomic proposition (AP) `:ap` to its labeling object.
  In a **concrete** problem, values are typically sets (e.g. LazySets / Dionysos sets) over the state space.
  In an **abstract** problem, values are typically collections of abstract states (e.g. `Vector{Int}`).

- `ap_semantics::Dict{Symbol, Any}`:
  Per-AP semantics used when converting set labels to abstract labels.
  Values: `Dionysos.Domain.INNER` or `Dionysos.Domain.OUTER`.

- `strict_spot::Bool`:
  If `true`, enforce Spot-style strict AP semantics / alphabet handling (useful to catch missing APs).
  If `false`, be permissive (e.g. treat missing APs as false / ignore absent labels depending on your pipeline).

This problem aims to synthesize a controller such that the generated trajectory satisfies the co-safe LTL
formula, i.e. it reaches an accepting condition in finite time.
"""
mutable struct CoSafeLTLProblem{S, XI, SPEC, LAB} <: ProblemType
    system::S
    initial_set::XI
    spec::SPEC

    # unified labeling container:
    labeling::Dict{Symbol, LAB}   # Symbol => LazySet (concrete) or Vector{Int} (abstract)

    ap_semantics::Dict{Symbol, Any}  # Symbol => DO.INNER / DO.OUTER
    strict_spot::Bool
end

struct Infinity <: Real end
Base.isfinite(::Infinity) = false

@recipe function f(problem::EmptyProblem; domain_color = :gray, region_color = :lightgray)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Region"
        color := region_color
        problem.region
    end
end

@recipe function f(
    problem::OptimalControlProblem;
    domain_color = :gray,
    initial_set_color = :green,
    target_set_color = :red,
)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
    @series begin
        label := "Target set"
        color := target_set_color
        problem.target_set
    end
end

@recipe function f(
    problem::SafetyProblem;
    domain_color = :gray,
    initial_set_color = :green,
    safe_set_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Safe set"
        color := safe_set_color
        problem.safe_set
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
end

@recipe function f(
    problem::CoSafeLTLProblem;
    domain_color = :gray,
    initial_set_color = :green,
    ap_colors = Dict{Symbol, Any}(),
    obs_color = :red,
)
    # --------------------
    # Domain
    # --------------------
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end

    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end

    for (ap, setX) in problem.labeling
        color_ap = haskey(ap_colors, ap) ? ap_colors[ap] : :blue
        @series begin
            label := String(ap)
            color := color_ap
            setX
        end
    end
end

export OptimalControlProblem
export SafetyProblem
export Infinity
