using StaticArrays, JuMP

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

function build_heuristic_controller(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    state_grid::DO.GridFree,
    input_grid::DO.GridFree,
    time_step::Real,
    jacobian_bound;
    approx_mode = AB.UniformGridAbstraction.GROWTH,
    sparse_input = true,
    print_level = 1,
    return_optimizer = false,
)
    # Step 1: Abstraction
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), time_step)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), approx_mode)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)

    MOI.optimize!(optimizer)

    # Step 2: Determinize symbolic model
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

    determinized_abstract_system = SY.determinize_symbolic_model(abstract_system)

    # Step 3: Solve new optimal control problem on determinized system
    new_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        new_optimizer,
        MOI.RawOptimizerAttribute("abstract_system"),
        determinized_abstract_system,
    )
    MOI.set(
        new_optimizer,
        MOI.RawOptimizerAttribute("discrete_time_system"),
        discrete_time_system,
    )
    MOI.set(new_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(new_optimizer, MOI.RawOptimizerAttribute("sparse_input"), sparse_input)
    MOI.set(new_optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)

    MOI.optimize!(new_optimizer)

    return return_optimizer ? new_optimizer :
           MOI.get(new_optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))
end

include("../problems/path_planning.jl")

concrete_problem = PathPlanning.problem(; simple = true)
state_grid = DO.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2))
input_grid = DO.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3))

value_fun = build_heuristic_controller(
    concrete_problem,
    state_grid,
    input_grid,
    0.3,
    PathPlanning.jacobian_bound();
    return_optimizer = false,
)

x0 = SVector(0.4, 0.4, 0.0)
println("Heuristic value for x₀ = $x0 is: ", value_fun(x0))
