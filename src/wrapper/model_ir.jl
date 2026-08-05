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
        MOI.FEASIBILITY_SENSE,
        nothing,
        0,
    )
end

function Base.empty!(ir::ModelIR)
    empty!(ir.variables)
    empty!(ir.dynamics)
    ir.time_domain = UNKNOWN
    empty!(ir.obstacles)
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
    if all(isnothing, ir.dynamics)
        error("Missing dynamics. i.e. ∂(x) = f(x, u) or Δ(x) = f(x, u)")
    end

    # Everything mentioned by some dynamics right-hand side.
    used = Set{Int}()
    for d in ir.dynamics
        d === nothing || _collect_variables!(used, d)
    end

    unused = Int[]
    state_index = 0
    input_index = 0
    for (i, v) in enumerate(ir.variables)
        if ir.dynamics[i] !== nothing
            v.role = STATE
            state_index += 1
            v.index = state_index
        else
            v.role = INPUT
            # A variable that is referenced nowhere and constrained nowhere cannot influence
            # the problem, so treating it as an input is never what the user meant.
            if !(i in used) && !_has_start(v) && !_has_target(v)
                push!(unused, i)
            end
            input_index += 1
            v.index = input_index
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

    return ir
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
