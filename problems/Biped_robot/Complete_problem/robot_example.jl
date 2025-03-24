using MathematicalSystems, StaticArrays, Plots
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
        "problems/Biped_robot/Complete_problem",
        "robot_problem.jl",
    ),
)

#######################################################
################### File Parameters ###################
#######################################################
filename_save = joinpath(@__DIR__, "Abstraction_solver.jld2")
do_empty_optim = false
verify_save = false

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

x0 = SVector{n_state, Float64}(ones(n_state).*(0)) # x0 in our case is a cell center !
hx = SVector{n_state, Float64}(fill(0.3, n_state)) # Intentional big discretization step (otherwise way too many values and infinite optimize!)
state_grid = DO.GridFree(x0, hx)

u0 = SVector{n_input, Float64}(zeros(n_input))
hu = SVector{n_input, Float64}(fill(3.0, n_input)) # Intentional big discretization step (otherwise way too many values and infinite optimize!)
input_grid = DO.GridFree(u0, hu)

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)   
MOI.set(optimizer, MOI.Silent(), true)  
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)


### Optimize
if(do_empty_optim)
    #######################################################
    ################# Empty Optimisations #################
    #######################################################
    MOI.optimize!(optimizer)

    # Save the abstraction solver
    my_abstraction_solver = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"))
    start_time = time()

    # Open the file in write mode and save the data
    jldopen(filename_save, "w") do file
        file["my_abstraction_solver"] = my_abstraction_solver
    end  # This block ensures the file is closed after writing

    end_time = time()
    save_time = end_time - start_time
    @info("Time elapsed to save : $save_time")

    #######################################################
    ######### Verify no info is lost in the save ##########
    #######################################################
    if(verify_save)
        # In the abtgsraction_solver, we have in the results : discrete_time_system, abstract_system and abstraction_construction_time_sec
        # We don't care about the last one. discrete_time_system is the function x[k+1] = f(x[k], u[k])
        # So we only need to compare the abstract_system
        jldopen(filename_save, "r") do file
            reloaded_solver = file["my_abstraction_solver"]
            # Perform the equality check within the read block to keep it scoped
            equal_test = RobotProblem.deep_equal(
                MOI.get(reloaded_solver, MOI.RawOptimizerAttribute("abstract_system")),
                MOI.get(my_abstraction_solver, MOI.RawOptimizerAttribute("abstract_system"))
            )

            @assert equal_test "Both abstract systems are not equal. See the print to see which field is the source."
        end
    end
else
    # Reload the result
    file = jldopen(filename_save, "r")
    reloaded_solver = file["my_abstraction_solver"]
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"), reloaded_solver)

    #######################################################
    ################# Problem definition ##################
    #######################################################
    _I_ = UT.HyperRectangle(x0, x0) # We force the system to start in the cell in which x_0 is

    t_low = SVector{n_state,Float64}([0.1, 0.1, -0.1, -0.1, -0.8, -0.8, -0.8, -0.8])
    t_high = SVector{n_state,Float64}([0.5, 0.5, -0.5, -0.5, 0.8, 0.8, 0.8, 0.8])
    _T_ = UT.HyperRectangle(t_low, t_high) # TODO
      
    concrete_problem = Dionysos.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        Dionysos.Problem.Infinity(),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        concrete_problem,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.optimize!(optimizer)
    abstract_problem_time =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
    println("Time to solve the abstract problem: $(abstract_problem_time)")
    
    abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
    uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))
    
    nstep = 30 # correspond to 3sec
    function reached(x)
        if x ∈ concrete_problem.target_set
            return true
        else
            return false
        end
    end

    control_trajectory = ST.get_closed_loop_trajectory(
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
    );
    
end



