module UniformEllipsoidAbstraction

import Dionysos
const UT = Dionysos.Utils
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic

import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA
import Polyhedra
import MathematicalSystems as MS

using ProgressMeter
using JuMP
import MathOptInterface as MOI
using JLD2

export Optimizer

include("empty_problem.jl")
include("optimal_control_problem.jl")

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    abstraction_solver::Union{Nothing, OptimizerEmptyProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, MS.AbstractSystem, MS.AbstractMap}
    solve_time_sec::T
    print_level::Int

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, 0.0, 1)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(opt::Optimizer) = opt.abstraction_solver === nothing

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    return model.print_level = value ? 0 : 1
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    param_symbol = Symbol(param.name)

    if model.abstraction_solver === nothing
        model.abstraction_solver = OptimizerEmptyProblem()
    end

    if param_symbol == :abstract_system
        MOI.set(
            model.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            value,
        )
        return
    end

    if param_symbol == :concrete_problem
        # Assign appropriate control solver
        if isa(value, Dionysos.Problem.EmptyProblem)
            model.abstraction_solver = OptimizerEmptyProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                value,
            )
            model.control_solver = nothing  # No control solver
        elseif isa(value, Dionysos.Problem.OptimalControlProblem)
            model.control_solver = OptimizerOptimalControlProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        else
            error("Unsupported problem type: $(typeof(value))")
        end

        # Instantiate an abstraction_solver if it has not already been created
        if model.abstraction_solver.empty_problem === nothing
            empty_problem = Dionysos.Problem.EmptyProblem(value.system, nothing)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                empty_problem,
            )
        end
        return
    end

    # If the attribute exists in the main optimizer, set it there
    if hasfield(typeof(model), param_symbol)
        return setproperty!(model, param_symbol, value)
    end

    # Try setting it in the sub-solvers
    for solver in (model.abstraction_solver, model.control_solver)
        if solver !== nothing && hasfield(typeof(solver), param_symbol)
            return setproperty!(solver, param_symbol, value)
        end
    end

    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    param_symbol = Symbol(param.name)

    # First, check if the attribute exists directly in the main optimizer
    if hasfield(typeof(model), param_symbol)
        return getproperty(model, param_symbol)
    end

    # If not found, try to get it from the abstraction solver
    if model.abstraction_solver !== nothing &&
       hasfield(typeof(model.abstraction_solver), param_symbol)
        return getproperty(model.abstraction_solver, param_symbol)
    end

    # If not found, try to get it from the control solver
    if model.control_solver !== nothing &&
       hasfield(typeof(model.control_solver), param_symbol)
        return getproperty(model.control_solver, param_symbol)
    end

    # If the attribute is not recognized, raise an error
    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

MOI.get(opt::Optimizer, ::MOI.SolveTimeSec) = opt.solve_time_sec

function is_abstraction_computed(opt::Optimizer)
    return opt.abstraction_solver !== nothing && opt.abstraction_solver.abstract_system !== nothing
end

function reset!(optimizer::Optimizer)
    optimizer.concrete_controller = nothing
    optimizer.solve_time_sec = 0.0
    if optimizer.control_solver !== nothing
        reset!(optimizer.control_solver)
    end
    if optimizer.abstraction_solver !== nothing
        reset!(optimizer.abstraction_solver)
    end
    return optimizer
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Ensure the concrete problem is defined
    if optimizer.abstraction_solver === nothing
        error("The concrete problem is not defined.")
    end

    # Compute abstraction if not already done
    if !is_abstraction_computed(optimizer)
        MOI.set(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        MOI.optimize!(optimizer.abstraction_solver)
    end

    # If there's a control solver, optimize it
    if optimizer.control_solver !== nothing
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        abstract_system = MOI.get(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
        )
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            abstract_system,
        )
        transitionCost = optimizer.abstraction_solver.transitionCost
        costfun = (q::Int, s::Int) -> get(transitionCost, (q, s), Inf)
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_transition_cost"),
            costfun,
        )

        MOI.optimize!(optimizer.control_solver)
        abstract_controller = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_controller"),
        )
        abstract_value_function = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_value_function"),
        )
        transitionCont = MOI.get(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("transitionCont"),
        )
        optimizer.concrete_controller = solve_concrete_problem(
            optimizer.abstraction_solver.abstract_system,
            abstract_controller,
            transitionCont,
            abstract_value_function
        )
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function solve_concrete_problem(
    abstract_system::SY.SymbolicModelList,
    abstract_controller::MS.AbstractMap,
    transitionCont::Dict,
    abstract_value_function::Function;
    handle_out_of_domain::Function = (x, abs_sys) -> nothing,
    randomize = false
)
    k_abs = abstract_controller.h
    is_defined_q = q -> (q ∈ abstract_controller.X)

    # choose a valid (from,to,cont) with minimal value
    function pick_best_transition(x)
        qs = SY.get_abstract_states(abstract_system, x)
        if isempty(qs)
            xnew = handle_out_of_domain(x, abstract_system)
            xnew === nothing && return nothing
            qs = SY.get_abstract_states(abstract_system, xnew)
            isempty(qs) && return nothing
            x = xnew
        end

        best_q = nothing
        best_to = nothing
        best_cont = nothing
        best_val = Inf

        for q in qs
            is_defined_q(q) || continue
            to_list = k_abs(q)
            to_list === nothing && return nothing
            to = isempty(to_list) ? nothing : (randomize ? rand(to_list) : first(to_list))
            to === nothing && return nothing
     
            key = (q, to)
            haskey(transitionCont, key) || continue

            v = abstract_value_function(q)
            if v < best_val
                best_val = v
                best_q = q
                best_to = to
                best_cont = transitionCont[key]
            end
        end

        if best_q === nothing
            return nothing
        end
        return (x = x, from = best_q, to = best_to, cont = best_cont, val = best_val)
    end

    # concrete controller x -> u
    f = function (x)
        tr = pick_best_transition(x)
        tr === nothing && return nothing
        return MS.apply(tr.cont, tr.x)
    end

    # domain predicate (controller defined at x)
    Xx = PredicateDomain(x -> pick_best_transition(x) !== nothing)

    nx = SY.get_state_dim(abstract_system)
    nu = SY.get_input_dim(abstract_system)
    return MS.ConstrainedBlackBoxMap(nx, nu, f, Xx)
end

# --------------------------------------- #
# --- JLD2 Abstraction Export/Import ---- #
# --------------------------------------- #

"""
export_abstraction_jld2(opt::UniformEllipsoidAbstraction.Optimizer, filename)

Stores:
- abstract_system
- minimal metadata to reload and understand it
"""
function export_abstraction_jld2(opt::Optimizer, filename::AbstractString)
    abs_opt = opt.abstraction_solver
    abs_opt === nothing && error("No abstraction_solver in optimizer.")
    abs_sys = abs_opt.abstract_system
    abs_sys === nothing && error("No abstract_system computed yet.")

    jldopen(filename, "w") do f
        f["format_version"] = 1
        f["abstract_system"] = abs_sys

        # store ellipsoid-abstraction-specific parameters that affect transitions
        f["params"] = (
            incl_mode = abs_opt.incl_mode,
            P = abs_opt.P,
            Pm = abs_opt.Pm,
            R = abs_opt.R,
        )
    end
    return nothing
end

"""
import_abstraction_jld2(filename; opt=nothing)

Loads abstract_system into the optimizer's abstraction_solver.
Returns the optimizer.
"""
function import_abstraction_jld2(filename::AbstractString; opt::Union{Nothing, Optimizer}=nothing)
    opt === nothing && (opt = MOI.instantiate(Optimizer))
    opt.abstraction_solver === nothing && (opt.abstraction_solver = OptimizerEmptyProblem())

    jldopen(filename, "r") do f
        v = f["format_version"]
        v == 1 || error("Unsupported abstraction file format_version=$v")

        abs_sys = f["abstract_system"]
        opt.abstraction_solver.abstract_system = abs_sys

        # optional: reload params if you want them accessible
        if haskey(f, "params")
            p = f["params"]
            # best-effort restore
            try opt.abstraction_solver.incl_mode = p.incl_mode catch end
            try opt.abstraction_solver.P = p.P catch end
            try opt.abstraction_solver.Pm = p.Pm catch end
            try opt.abstraction_solver.R = p.R catch end
        end
    end

    return opt
end

end # module