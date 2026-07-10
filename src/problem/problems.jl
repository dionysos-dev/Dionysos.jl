struct Infinity end
Base.isfinite(::Infinity) = false

discretize_time(time::Real, Δt::Real; round_up = true) =
    round_up ? ceil(Int, time / Δt) : floor(Int, time / Δt)
discretize_time(::Infinity, Δt) = Infinity()

"""
    ProblemType

An abstract type that represents a generic control problem.  
All concrete problem types (e.g., [`AlternatingSimulationProblem{S,X}`](@ref), [`OptimalControlProblem`](@ref), [`SafetyProblem`](@ref), [`CoSafeLTLProblem`](@ref)) should subtype `ProblemType`.
"""
abstract type ProblemType end

"""
    AlternatingSimulationProblem{S, X} <: ProblemType

A problem type used to construct a **sound abstraction** of a dynamical system.

- `S`: The system to abstract (continuous or discrete-time).
- `X`: The region of interest (e.g., a subset of the state space).

This problem encodes no control objective. It is intended for generating symbolic models that can later be reused by other solvers.
"""
struct AlternatingSimulationProblem{S, X} <: ProblemType
    system::S
    region::X
end

"""
    BisimulationQuotientProblem{S, X, R} <: ProblemType

A problem type used to construct a finite bisimulation (exact equivalence abstraction) quotient induced.

# Fields
- `system`: the switched system to abstract.
- `region`: the state-space region of interest `X`.
- `observation_regions`: the regions of interest used to define the observation map.

This problem encodes no control objective. It is intended for generating symbolic models that can later be reused by other solvers.
"""
struct BisimulationQuotientProblem{S, X, R} <: ProblemType
    system::S
    region::X
    observation_regions::R
end

"""
    OptimalControlProblem{S, XI, XT, XC, TC, T} <: ProblemType

Encodes a **reach-avoid optimal control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached.
- `XC`: A state cost function or structure.
- `TC`: A transition cost function or structure.
- `T`: Satisfy the property in at most time T.

This problem aims to find a control strategy that reaches the target set from the initial set, minimizing the accumulated cost over time.
"""
struct OptimalControlProblem{S, XI, XT, XC, TC, T} <: ProblemType
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
end

OptimalControlProblem(system, initial_set, target_set, state_cost, transition_cost) =
    OptimalControlProblem(
        system,
        initial_set,
        target_set,
        state_cost,
        transition_cost,
        Infinity(),
    )

"""
    SafetyProblem{S, XI, XS, T} <: ProblemType

Encodes a **safety control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XS`: The safe set in which the system must remain.
- `T`: Satisfy the property for at least time T.

This problem aims to synthesize a controller that ensures the system remains within the safe set for the entire duration of the time horizon.
"""
struct SafetyProblem{S, XI, XS, T} <: ProblemType
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

SafetyProblem(system, initial_set, safe_set) =
    SafetyProblem(system, initial_set, safe_set, Infinity())

"""
    ReachAndStayProblem{S, XI, XT, XS, T} <: ProblemType

Encodes a **reach-and-stay control problem** (eventually always).

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached and stayed in.
- `XS`: The safe set in which the system must remain during the approach.
- `T`: Satisfy the property for at least time T.

This problem aims to synthesize a controller that drives the system from the initial set
into the target set and keeps it there indefinitely, while remaining within the safe set
during the approach phase.
"""
struct ReachAndStayProblem{S, XI, XT, XS, T} <: ProblemType
    system::S
    initial_set::XI
    target_set::XT
    safe_set::XS
    time::T
end

ReachAndStayProblem(system, initial_set, target_set, safe_set) =
    ReachAndStayProblem(system, initial_set, target_set, safe_set, Infinity())

"""
    CoSafeLTLProblem{S, XI, SPEC, LAB} <: ProblemType

Encodes a **co-safe LTL control problem**.

- `S`: The system to control.
- `XI`: The initial set of states.
- `SPEC`: The co-safe LTL specification object (an automaton/monitor wrapper).
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

- `ap_semantics::Dict{Symbol, UT.INCL_MODE}`:
  Per-AP semantics used when converting set labels to abstract labels
  (`UT.INNER` or `UT.OUTER`; also reachable as `Mapping.INNER`/`Mapping.OUTER`).


This problem aims to synthesize a controller such that the generated trajectory satisfies the co-safe LTL
formula, i.e. it reaches an accepting condition in finite time.
"""
struct CoSafeLTLProblem{S, XI, SPEC, LAB} <: ProblemType
    system::S
    initial_set::XI
    spec::SPEC

    # unified labeling container:
    labeling::Dict{Symbol, LAB}   # Symbol => LazySet (concrete) or Vector{Int} (abstract)

    ap_semantics::Dict{Symbol, UT.INCL_MODE}
end

# Parametric structs do not convert non-parametric fields, so accept any
# INCL_MODE-valued dict (e.g. the common `Dict{Symbol, Any}`) and convert.
function CoSafeLTLProblem(
    system,
    initial_set,
    spec,
    labeling::Dict{Symbol, LAB},
    ap_semantics::AbstractDict{Symbol},
) where {LAB}
    return CoSafeLTLProblem(
        system,
        initial_set,
        spec,
        labeling,
        convert(Dict{Symbol, UT.INCL_MODE}, ap_semantics),
    )
end

trajectory_success(problem::ProblemType, traj::ST.Trajectory) = false

function trajectory_success(problem::OptimalControlProblem, traj::ST.Trajectory)
    isempty(traj.seq) && return false

    return first(traj.seq) ∈ problem.initial_set &&
           any(x -> x ∈ problem.target_set, traj.seq)
end

function trajectory_success(problem::SafetyProblem, traj::ST.Trajectory)
    isempty(traj.seq) && return false

    return first(traj.seq) ∈ problem.initial_set && all(x -> x ∈ problem.safe_set, traj.seq)
end

function trajectory_success(problem::ReachAndStayProblem, traj::ST.Trajectory)
    xs = traj.seq
    first(xs) ∈ problem.initial_set || return false
    all(x -> x ∈ problem.safe_set, xs) || return false

    return any(k -> all(x -> x ∈ problem.target_set, xs[k:end]), eachindex(xs))
end

function trajectory_success(problem::CoSafeLTLProblem, traj::ST.Trajectory)
    isempty(traj.seq) && return false

    # Placeholder until monitor/spec trajectory evaluation is implemented.
    return false
end

function discretize_problem(
    problem::ProblemType,
    Δt::Float64;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    return error("implement `discretize_problem` for $(typeof(problem))")
end

function discretize_problem(
    problem::OptimalControlProblem,
    Δt::Float64;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    discrete_system =
        ST.discretize_continuous_system(problem.system, Δt; num_substeps = num_substeps)

    return OptimalControlProblem(
        discrete_system,
        problem.initial_set,
        problem.target_set,
        problem.state_cost,
        problem.transition_cost,
        discretize_time(problem.time, Δt; round_up = false),
    )
end

function discretize_problem(
    problem::SafetyProblem,
    Δt::Float64;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    discrete_system =
        ST.discretize_continuous_system(problem.system, Δt; num_substeps = num_substeps)

    return SafetyProblem(
        discrete_system,
        problem.initial_set,
        problem.safe_set,
        discretize_time(problem.time, Δt; round_up = true),
    )
end

function discretize_problem(
    problem::ReachAndStayProblem,
    Δt::Float64;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    discrete_system =
        ST.discretize_continuous_system(problem.system, Δt; num_substeps = num_substeps)

    return ReachAndStayProblem(
        discrete_system,
        problem.initial_set,
        problem.target_set,
        problem.safe_set,
        discretize_time(problem.time, Δt; round_up = true),
    )
end

function discretize_problem(
    problem::CoSafeLTLProblem,
    Δt::Float64;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    discrete_system =
        ST.discretize_continuous_system(problem.system, Δt; num_substeps = num_substeps)

    return CoSafeLTLProblem(
        discrete_system,
        problem.initial_set,
        problem.spec,
        problem.labeling,
        problem.ap_semantics,
    )
end

@recipe function f(
    problem::AlternatingSimulationProblem;
    domain_color = :gray,
    region_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end
    @series begin
        label := "Region"
        color := region_color
        problem.region
    end
end

@recipe function f(
    problem::BisimulationQuotientProblem;
    region_color = :gray,
    region_alpha = 0.15,
    observation_region_alpha = 0.4,
    observation_colors = [:red, :green, :orange, :purple, :brown],
    plot_region = true,
)
    if plot_region
        @series begin
            label := "Region"
            color := region_color
            fillalpha := region_alpha
            problem.region
        end
    end
    for (i, R) in enumerate(problem.observation_regions)
        @series begin
            label := "O $i"
            color := observation_colors[mod1(i, length(observation_colors))]
            fillalpha := observation_region_alpha
            R
        end
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
        MS.stateset(problem.system)
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
        MS.stateset(problem.system)
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
    problem::ReachAndStayProblem;
    domain_color = :gray,
    safe_set_color = :lightgray,
    target_set_color = :red,
    initial_set_color = :green,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end

    @series begin
        label := "Safe set"
        color := safe_set_color
        problem.safe_set
    end

    @series begin
        label := "Target set"
        color := target_set_color
        problem.target_set
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
        MS.stateset(problem.system)
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
export ReachAndStayProblem
export CoSafeLTLProblem
export Infinity
