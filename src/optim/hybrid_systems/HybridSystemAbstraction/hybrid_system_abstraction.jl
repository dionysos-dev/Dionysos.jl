export HybridSystemAbstraction

module HybridSystemAbstraction
using MathOptInterface
using HybridSystems
import MathematicalSystems
MS = MathematicalSystems
using StaticArrays: SVector
using Dionysos

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems

export Optimizer

include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    abstraction_solver::Union{Nothing, OptimizerAlternatingSimulationProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, ST.AbstractContinuousController}
    solve_time_sec::T
    print_level::Int

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, 0.0, 1)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.abstraction_solver === nothing

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    return model.print_level = value ? 0 : 1
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    param_symbol = Symbol(param.name)

    if model.abstraction_solver === nothing
        model.abstraction_solver = OptimizerAlternatingSimulationProblem()
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
        if isa(value, PR.AlternatingSimulationProblem)
            model.abstraction_solver = OptimizerAlternatingSimulationProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                value,
            )
            model.control_solver = nothing  # No control solver
        elseif isa(value, PR.OptimalControlProblem)
            model.control_solver = OptimizerOptimalControlProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        elseif isa(value, PR.SafetyProblem)
            model.control_solver = OptimizerSafetyProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        elseif isa(value, PR.CoSafeLTLProblem)
            model.control_solver = OptimizerCoSafeLTLProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        else
            error("Unsupported problem type: $(typeof(value))")
        end

        # Instantiate an abstraction_solver if it has not already been created
        if model.abstraction_solver.alternating_simulation_problem === nothing
            alternating_simulation_problem =
                PR.AlternatingSimulationProblem(value.system, nothing)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                alternating_simulation_problem,
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

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function is_abstraction_computed(optimizer::Optimizer)
    return optimizer.abstraction_solver !== nothing &&
           optimizer.abstraction_solver.abstract_system !== nothing
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
        MOI.optimize!(optimizer.control_solver)
        abstract_controller = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_controller"),
        )
        optimizer.concrete_controller =
            HybridQuantizedStaticController(abstract_system, abstract_controller)
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

# ================================================================
# Concrete controller synthesis
# ================================================================

struct HybridQuantizedStaticController{AS, AC} <: ST.AbstractContinuousController
    abstract_system::AS
    abstract_controller::AC
end
ST.controller_kind(::HybridQuantizedStaticController) = ST.StaticKind()
ST.domain(ctrl::HybridQuantizedStaticController) = ctrl.abstract_system
ST.initial_state(::HybridQuantizedStaticController) = nothing
ST.update_state(::HybridQuantizedStaticController, q, y) = nothing

function ST.is_defined(ctrl::HybridQuantizedStaticController, q, aug_state)
    abs_q = SY.get_abstract_state(ctrl.abstract_system, aug_state)
    abs_q === nothing && return false
    return ST.is_defined(ctrl.abstract_controller, nothing, abs_q)
end

function ST.output_control(ctrl::HybridQuantizedStaticController, q, aug_state)
    abs_q = SY.get_abstract_state(ctrl.abstract_system, aug_state)
    abs_q === nothing && return nothing

    u_abs = ST.output_control(ctrl.abstract_controller, nothing, abs_q)
    u_abs === nothing && return nothing

    (_, _, k) = aug_state
    if SY.is_switching_input(ctrl.abstract_system.input_mapping, u_abs)
        transition_id = ctrl.abstract_system.input_mapping.global_to_switching[u_abs]
        return ctrl.abstract_system.input_mapping.switch_labels[transition_id]
    else
        return SY.get_concrete_input(ctrl.abstract_system, u_abs, k)
    end
end

# ================================================================
# Closed-loop simulation utilities
# ================================================================

function get_next_aug_state(hs::HybridSystem, aug_state, u, time_is_active, tstep, map_sys)
    (x, t, k) = aug_state

    if isa(u, AbstractString) && occursin("SWITCH", u)
        m = match(r"SWITCH (\d+) -> (\d+)", String(u))
        if m === nothing
            error("Unrecognized SWITCH label format: $u")
        end
        source_mode = parse(Int, m.captures[1])
        target_mode = parse(Int, m.captures[2])

        transitions = collect(HybridSystems.transitions(hs.automaton))
        transition_id = findfirst(
            tr ->
                HybridSystems.source(hs.automaton, tr) == source_mode &&
                HybridSystems.target(hs.automaton, tr) == target_mode,
            transitions,
        )
        if transition_id === nothing
            error("Transition not found for $u")
        end
        transition = transitions[transition_id]

        reset_map = HybridSystems.resetmap(hs, transition)
        augmented_source_state = vcat(x, t)
        reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
        next_x = reset_result[1:(end - 1)]
        next_t = reset_result[end]
        next_k = target_mode
        next_t = round(next_t; digits = 10)
        return (next_x, next_t, next_k)
    else
        next_t = time_is_active ? t + tstep : 0.0
        next_t = round(next_t; digits = 10)
        next_x = map_sys(x, u, tstep)
        return (next_x, next_t, k)
    end
end

function get_closed_loop_trajectory(
    hs::HybridSystem,
    controller::ST.AbstractController,
    tsteps,
    aug_state_0,
    nstep;
    stopping = x -> false,
)
    q = ST.initial_state(controller)
    aug_state_traj, u_traj = [aug_state_0], Any[]
    aug_state = aug_state_0

    nmodes = HybridSystems.nmodes(hs.automaton)
    dynamics = [HybridSystems.mode(hs, k).systems[1].f for k in 1:nmodes]
    maps_sys = [ST.simulate_control_map(dynamics[i]) for i in 1:nmodes]
    times_is_active =
        [([1.0;;] == HybridSystems.mode(hs, k).systems[2].A) for k in 1:nmodes]

    for _ in 1:nstep
        stopping(aug_state) && break

        u = ST.output_control(controller, q, aug_state)
        u === nothing && break

        (_, _, k) = aug_state
        next_aug_state =
            get_next_aug_state(hs, aug_state, u, times_is_active[k], tsteps[k], maps_sys[k])

        q = ST.update_state(controller, q, aug_state)
        aug_state = next_aug_state

        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end

    return aug_state_traj, u_traj
end

end # module
