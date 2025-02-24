module TestMain
using Test     #src

using StaticArrays, Plots

# At this point, we import Dionysos.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# And the file defining the hybrid system for this problem
include("../../problems/path_planning.jl")

# ### Definition of the problem

# Now we instantiate the problem using the function provided by [PathPlanning.jl](@__REPO_ROOT_URL__/problems/PathPlanning.jl) 
concrete_problem = PathPlanning.problem(; simple = true)
concrete_system = concrete_problem.system

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
x0 = SVector(0.0, 0.0, 0.0)
h = SVector(0.2, 0.2, 0.2)
state_grid = DO.GridFree(x0, h)

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0, 0.0)
h = SVector(0.3, 0.3)
input_grid = DO.GridFree(u0, h)

# We now solve the optimal control problem with the `Abstraction.UniformGridAbstraction.Optimizer`.

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    PathPlanning.jacobian_bound(),
)
MOI.optimize!(optimizer)

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

@testset "UniformGridAbstraction reachability" begin
    @test length(abstract_controller.data) == 19400 #src
end

# ### Trajectory display
# We choose a stopping criterion `reached` and the maximal number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space `_I_` defined previously.
nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end
x0 = SVector(0.4, 0.4, 0.0)
control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

no_plot = false
@static if get(ENV, "CI", "false") == "false" &&
           (isdefined(@__MODULE__, :no_plot) && no_plot == false)

    # Here we display the coordinate projection on the two first components of the state space along the trajectory.
    fig = plot(; aspect_ratio = :equal)
    # We display the concrete domain
    plot!(concrete_system.X; color = :yellow, opacity = 0.5)

    # We display the abstract domain
    plot!(abstract_system.Xdom; color = :blue, opacity = 0.5)

    # We display the concrete specifications
    plot!(concrete_problem.initial_set; color = :green, opacity = 0.2)
    plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2)

    # We display the abstract specifications
    plot!(
        SY.get_domain_from_states(abstract_system, abstract_problem.initial_set);
        color = :green,
    )
    plot!(
        SY.get_domain_from_states(abstract_system, abstract_problem.target_set);
        color = :red,
    )

    # We display the concrete trajectory
    plot!(control_trajectory; ms = 0.5)
end
end
# ### References
# 1. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
# 2. K. J. Aström and R. M. Murray, Feedback systems. Princeton University Press, Princeton, NJ, 2008.
