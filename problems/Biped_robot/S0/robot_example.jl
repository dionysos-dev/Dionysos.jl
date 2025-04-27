using MathematicalSystems, StaticArrays, Plots, LinearAlgebra
using JuMP
using JLD2

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems/Biped_robot/S0",
        "robot_problem.jl",
    ),
)

#######################################################
################### File Parameters ###################
#######################################################
First_part = false
Second_part = true

#######################################################
################### Optim Parameters ##################
#######################################################
concrete_problem = RobotProblem.problem(; tstep = 1e-1)
concrete_system = concrete_problem.system

### Set the optimizer
n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)
state_space = MathematicalSystems.stateset(concrete_system)
input_space = MathematicalSystems.inputset(concrete_system)
println("n_state: ", n_state)
println("n_input: ", n_input)
x0 = SVector{n_state, Float64}(zeros(n_state))

if(First_part)
    println("First part : ")
    println()
    
    filename_controller_1 = joinpath(@__DIR__, "Controller_step_1.jld2")
    file = jldopen(filename_controller_1, "r")
    reloaded_optimizer_1 = file["optimizer"]

    concrete_problem = MOI.get(reloaded_optimizer_1,  MOI.RawOptimizerAttribute("concrete_problem"))
    
    nstep = 300 # correspond to 30 sec
    function reached(x)
        if x ∈ concrete_problem.target_set
            return true
        else
            return false
        end
    end
    println(concrete_problem.target_set)

    #######################################################
    ################# Concrete Trajectory #################
    #######################################################
    control_trajectory = ST.get_closed_loop_trajectory(
        concrete_system,
        MOI.get(reloaded_optimizer_1, MOI.RawOptimizerAttribute("concrete_controller")),
        x0,
        nstep;
        stopping = reached,
    );
    println(control_trajectory)
    println()

end
if(Second_part)
    println("Second Part")
    println()

    #######################################################
    ################# Problem definition ##################
    #######################################################
    p0 = SVector{n_state,Float64}([-0.18042316250653628, 0.13622830985740686, 0.206712192483932, 0.0, 0.0, 0.0])

    filename_controller_2 = joinpath(@__DIR__, "Controller_step_2.jld2")
    file2 = jldopen(filename_controller_2, "r")
    reloaded_optimizer_2 = file2["optimizer"]
    reloaded_optimizer_2.handle_out_of_domain = 0

    concrete_problem = MOI.get(reloaded_optimizer_2,  MOI.RawOptimizerAttribute("concrete_problem"))
    
    nstep = 300 # correspond to 30 sec
    
    function reached(x)
        if x ∈ concrete_problem.target_set
            return true
        else
            return false
        end
    end

    # In our case, the state_space has been defined too large when we did the abstraction
    # We have to re-create the exact one
    half_step = [fill(π/180, 3)..., fill(0.075, 3)...]
    state_lower_bounds = [-12*π/180, 0, 0, -0.6, -0.15, -0.15] .- half_step
    state_upper_bounds = [0, 12*π/180, 14*π/180, 0.15, 0.6, 0.6] .+ half_step

    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    concrete_control_trajectory = ST.get_closed_loop_trajectory(
        concrete_system,
        MOI.get(reloaded_optimizer_2, MOI.RawOptimizerAttribute("concrete_controller")),
        p0,
        nstep;
        stopping = reached,
        handle_out_of_domain = reloaded_optimizer_2.handle_out_of_domain,
        state_space = state_space
    );

    println(concrete_control_trajectory)
    println()
    
end
