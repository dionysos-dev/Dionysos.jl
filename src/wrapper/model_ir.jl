# ----------------------------------------------------------------------------------------
# ModelIR — the parsed JuMP model.
#
# This is the seam of the whole front-end: `parse.jl` fills it, `lower_*.jl` consume it, and
# nothing in it depends on Symbolics. Validation runs once here, on a complete model, instead
# of incrementally inside `MOI.add_constraint`.
# ----------------------------------------------------------------------------------------

"""
    TimeDomain

Whether the model was written with continuous-time (`∂`, `CONTINUOUS`) or discrete-time
(`Δ`, `DISCRETE`) dynamics. `UNKNOWN` until the first dynamics constraint fixes it; mixing the
two is an error.
"""
@enum TimeDomain UNKNOWN CONTINUOUS DISCRETE

# The dynamics of one scope: state index → the expression driving it.
const DynamicsMap = Dict{Int, MOI.ScalarNonlinearFunction}

"""
    ModeIR

One discrete mode: its own dynamics, its own bound overrides (a mode may restrict the state or
input set — the thermostat's heater is off in one mode and throttled in the other), its own
specifications, and its own solver options.
"""
mutable struct ModeIR
    id::Int
    dynamics::DynamicsMap
    user_dynamics::Any
    lower::Dict{Int, Float64}
    upper::Dict{Int, Float64}
    specs::Vector{SpecEntry}
    attributes::Vector{Pair{String, Any}}
end

ModeIR(id::Int) = ModeIR(
    id,
    DynamicsMap(),
    nothing,
    Dict{Int, Float64}(),
    Dict{Int, Float64}(),
    SpecEntry[],
    Pair{String, Any}[],
)

"""
    TransitionIR

One switch: the guard bounds that enable it and the reset map applied when it is taken.
"""
mutable struct TransitionIR
    id::Int
    source::Int
    target::Int
    switching::Any
    guard_lower::Dict{Int, Float64}
    guard_upper::Dict{Int, Float64}
    resets::DynamicsMap
end

TransitionIR(id::Int, source::Int, target::Int, switching) = TransitionIR(
    id,
    source,
    target,
    switching,
    Dict{Int, Float64}(),
    Dict{Int, Float64}(),
    DynamicsMap(),
)

"""
    ModelIR

Dependency-free description of a parsed JuMP model: variables and their inferred roles, the
per-variable dynamics expressions, obstacles, and the objective.

Fields are populated by `parse.jl` as JuMP hands constraints over, then validated in one pass
by [`infer_roles!`](@ref) before lowering.
"""
mutable struct ModelIR
    variables::Vector{VariableInfo}
    dynamics::Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}
    time_domain::TimeDomain
    obstacles::Vector{Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}}
    specs::Vector{SpecEntry}
    labels::Vector{LabelEntry}
    specification::Any
    modes::Dict{Int, ModeIR}
    transitions::Dict{Int, TransitionIR}
    # A Julia function supplied instead of written equations (issue #510).
    user_dynamics::Any
    horizon::Union{Nothing, Float64}
    objective_sense::MOI.OptimizationSense
    objective_function::Union{Nothing, MOI.AbstractScalarFunction}
    nonlinear_index::Int
end

function ModelIR()
    return ModelIR(
        VariableInfo[],
        Union{Nothing, MOI.ScalarNonlinearFunction}[],
        UNKNOWN,
        Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}[],
        SpecEntry[],
        LabelEntry[],
        nothing,
        Dict{Int, ModeIR}(),
        Dict{Int, TransitionIR}(),
        nothing,
        nothing,
        MOI.FEASIBILITY_SENSE,
        nothing,
        0,
    )
end

"""
    is_hybrid(ir) -> Bool

Whether the model declared any mode.
"""
is_hybrid(ir::ModelIR) = !isempty(ir.modes)

"""
    has_user_dynamics(ir) -> Bool

Whether any dynamics were supplied as a Julia function rather than written as equations.
"""
has_user_dynamics(ir::ModelIR) =
    ir.user_dynamics !== nothing || any(m -> m.user_dynamics !== nothing, values(ir.modes))

"""
    mode_ids(ir) -> Vector{Int}

Mode identifiers in declaration order. They index the `HybridSystems` automaton directly, so
they must be `1:n`; [`infer_roles!`](@ref) checks that.
"""
mode_ids(ir::ModelIR) = sort!(collect(keys(ir.modes)))

# The mode / transition record for `id`, created on first mention. Scopes are discovered from
# the constraints that carry them, so there is no separate declaration channel.
mode!(ir::ModelIR, id::Int) = get!(() -> ModeIR(id), ir.modes, id)

function transition!(ir::ModelIR, scope)
    return get!(ir.transitions, scope.id) do
        return TransitionIR(scope.id, scope.source, scope.target, scope.switching)
    end
end

function Base.empty!(ir::ModelIR)
    empty!(ir.variables)
    empty!(ir.dynamics)
    ir.time_domain = UNKNOWN
    empty!(ir.obstacles)
    empty!(ir.specs)
    empty!(ir.labels)
    ir.specification = nothing
    empty!(ir.modes)
    empty!(ir.transitions)
    ir.user_dynamics = nothing
    ir.horizon = nothing
    ir.objective_sense = MOI.FEASIBILITY_SENSE
    ir.objective_function = nothing
    ir.nonlinear_index = 0
    return ir
end

Base.isempty(ir::ModelIR) = isempty(ir.variables)

function add_variable!(ir::ModelIR)
    push!(ir.variables, VariableInfo())
    push!(ir.dynamics, nothing)
    return MOI.VariableIndex(length(ir.variables))
end

# Fresh index for the constraints the wrapper consumes structurally (dynamics, start, final):
# they are recorded in the IR rather than stored as constraints, so they only need a unique
# handle to hand back to JuMP.
function next_constraint_index!(ir::ModelIR, func, set)
    ir.nonlinear_index += 1
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(ir.nonlinear_index)
end

"""
    infer_roles!(ir::ModelIR)

Assign a [`VariableRole`](@ref) to every variable and validate the model as a whole.

Rules (`src/wrapper/README.md` §3): a variable carrying `∂`/`Δ` dynamics is a `STATE` (I1);
any other variable appearing on a dynamics right-hand side is an `INPUT` (I3); a variable that
appears nowhere is an error (I4), because it would otherwise be silently discretized as an
input and enlarge the abstraction.
"""
function infer_roles!(ir::ModelIR)
    _validate_modes(ir)

    # A variable is a state if it carries dynamics anywhere: at the top level, or in *any*
    # mode. A hybrid model may well leave a state undriven in one mode.
    has_dynamics = falses(length(ir.variables))
    for (i, d) in enumerate(ir.dynamics)
        d === nothing || (has_dynamics[i] = true)
    end
    for m in values(ir.modes), i in keys(m.dynamics)
        has_dynamics[i] = true
    end

    # Dynamics supplied as a Julia function leave no expressions to infer from, so the states
    # are named with `set_role!` instead and the checks below stand down.
    supplied = has_user_dynamics(ir)

    if !any(has_dynamics) && !supplied
        error("Missing dynamics. i.e. ∂(x) = f(x, u) or Δ(x) = f(x, u)")
    end

    # Everything mentioned by some dynamics or reset right-hand side.
    used = Set{Int}()
    for d in ir.dynamics
        d === nothing || _collect_variables!(used, d)
    end
    for m in values(ir.modes), d in values(m.dynamics)
        _collect_variables!(used, d)
    end
    for t in values(ir.transitions), d in values(t.resets)
        _collect_variables!(used, d)
    end
    for t in values(ir.transitions)
        union!(used, keys(t.guard_lower))
        union!(used, keys(t.guard_upper))
    end

    unused = Int[]
    state_index = 0
    input_index = 0
    for (i, v) in enumerate(ir.variables)
        declared = v.declared_role
        v.role = declared !== nothing ? declared : (has_dynamics[i] ? STATE : INPUT)
        if v.role === INPUT
            # A variable that is referenced nowhere and constrained nowhere cannot influence
            # the problem, so treating it as an input is never what the user meant. With
            # supplied dynamics there is nothing to be referenced by, so the check cannot run.
            if declared === nothing &&
               !supplied &&
               !(i in used) &&
               !_has_start(v) &&
               !_has_target(v)
                push!(unused, i)
            end
            input_index += 1
            v.index = input_index
        else
            state_index += 1
            v.index = state_index
        end
    end

    if !isempty(unused)
        names = join((describe(ir.variables[i], i) for i in unused), ", ")
        error(
            "Unused variable(s): $names. They appear in no dynamics and carry no `start`/`final` " *
            "constraint, so they would be silently discretized as inputs and enlarge the " *
            "abstraction. Remove them, use them in the dynamics, or set their role explicitly.",
        )
    end

    # A clock is a state, but not a coordinate of the physical `x`, so it is re-labelled after
    # the state/input split and drops out of `state_indices`.
    detect_clock!(ir)

    _validate_bounds(ir)
    _validate_objective(ir)
    return ir
end

# Mode ids index the `HybridSystems` automaton directly, so they must be a contiguous 1:n.
# `add_mode!` allocates them that way; a gap means modes were declared on two different JuMP
# models, or one was created and never used.
function _validate_modes(ir::ModelIR)
    isempty(ir.modes) && return
    ids = mode_ids(ir)
    ids == collect(1:length(ids)) || error(
        "Mode ids must be contiguous from 1, got $ids. Declare every mode on the same model " *
        "with `@mode` / `add_mode!`, and give each one at least one constraint.",
    )
    for t in values(ir.transitions)
        (t.source in ids && t.target in ids) || error(
            "Transition $(t.id) runs between modes $(t.source) → $(t.target), but the model " *
            "only has modes $ids.",
        )
    end
    return
end

# The abstraction discretizes X and U, so every state and input must be boxed. Without this
# the unbounded coordinate reaches `UT.box` as ±Inf and throws a NaN-radius assertion from
# deep inside LazySets, saying nothing about which variable is at fault.
function _validate_bounds(ir::ModelIR)
    unbounded = findall(eachindex(ir.variables)) do i
        v = ir.variables[i]
        return !(isfinite(v.lower) && isfinite(v.upper))
    end
    isempty(unbounded) && return
    names = join((describe(ir.variables[i], i) for i in unbounded), ", ")
    return error(
        "Unbounded variable(s): $names. Every state and input must have finite bounds, " *
        "because the abstraction discretizes them; declare them as " *
        "`@variable(model, lo <= x <= hi)`.",
    )
end

function _validate_objective(ir::ModelIR)
    ir.objective_sense === MOI.FEASIBILITY_SENSE && return
    return error(
        "`@objective` is not supported yet: the wrapper would silently ignore it. A control " *
        "cost is set with `set_attribute(model, \"transition_cost\", (x, u) -> …)`; the " *
        "specification itself is expressed with constraints (`Final`, `Always`, …).",
    )
end

"""
    state_indices(ir) -> Vector{Int}

MOI indices of the state variables, in declaration order — this is the order of `x` in the
lowered system.
"""
state_indices(ir::ModelIR) =
    findall(i -> ir.variables[i].role === STATE, eachindex(ir.variables))

"""
    input_indices(ir) -> Vector{Int}

MOI indices of the input variables, in declaration order.
"""
input_indices(ir::ModelIR) =
    findall(i -> ir.variables[i].role === INPUT, eachindex(ir.variables))

function Base.show(io::IO, ir::ModelIR)
    nstate = count(!isnothing, ir.dynamics)
    print(
        io,
        "ModelIR($(length(ir.variables)) variables, $nstate with dynamics, ",
        "$(ir.time_domain), $(length(ir.obstacles)) obstacle(s))",
    )
    return
end
