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
        "problems/Biped_robot/Simplified_problem",
        "robot_problem.jl",
    ),
)

function get_abstract_closed_loop_trajectory(abstract_system, abstract_controller, source, nstep; stopping = (s)->false)
    state_traj, input_traj = [source], []

    for _ in 1:nstep
        stopping(source) && break
        abstract_inputs = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(abstract_inputs)
            @warn("Uncontrollable state: $source")
            return nothing
        end
        input = first(abstract_inputs)[1]
        targets = []
        UT.fix_and_eliminate_tail!(targets, abstract_system.autom.transitions, (source, input))
        source = first(targets)
        push!(state_traj, source)
        push!(input_traj, input)
    end
    return ST.Control_trajectory(ST.Trajectory(state_traj), ST.Trajectory(input_traj))
end

function get_concrete_trajectory(abstract_system, abstract_trajectory::ST.Control_trajectory)
    concrete_state_traj, concrete_input_traj = [], []
    for k in 1:ST.length(abstract_trajectory)
        abstarct_state = ST.get_state(abstract_trajectory, k)
        concrete_state = SY.get_concrete_state(abstract_system, abstarct_state)
        push!(concrete_state_traj, concrete_state)

        if k<ST.length(abstract_trajectory)
            abstract_input = ST.get_input(abstract_trajectory, k)
            concrete_input = SY.get_concrete_input(abstract_system, abstract_input)
            push!(concrete_input_traj, concrete_input)
        end
    end
    return ST.Control_trajectory(ST.Trajectory(concrete_state_traj), ST.Trajectory(concrete_input_traj))
end

#######################################################
################### File Parameters ###################
#######################################################
filename_save = joinpath(@__DIR__, "Abstraction_solver.jld2")
do_empty_optim = false
verify_save = false
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
hx = SVector{n_state, Float64}([fill(2*π/180, 3)..., fill(0.15, 3)...])
state_grid = DO.GridFree(x0, hx)

u0 = SVector{n_input, Float64}(zeros(n_input))
hu = SVector{n_input, Float64}(fill(1.0, n_input))
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
else
    if(First_part)
        println("First part : ")
        println()
        # Reload the result
        file = jldopen(filename_save, "r")
        reloaded_solver = file["my_abstraction_solver"]
        MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"), reloaded_solver)

        #######################################################
        ################# Problem definition ##################
        #######################################################
        _I_ = UT.HyperRectangle(x0, x0) # We force the system to start in the cell in which x_0 is

        t_low = SVector{n_state,Float64}([-12*π/180.0, 7*π/180.0, 8*π/180.0, -0.75, -0.3, -0.3])
        t_high = SVector{n_state,Float64}([-8*π/180.0, 9*π/180.0, 12*π/180.0, 0.3, 0.75, 0.75])
        _T_ = UT.HyperRectangle(t_low, t_high)
        
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
        abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
        concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
        controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
        uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))
        
        nstep = 300 # correspond to 30 sec
        function reached(x)
            if x ∈ concrete_problem.target_set
                return true
            else
                return false
            end
        end

        #######################################################
        ################# Concrete Trajectory #################
        #######################################################
        control_trajectory = ST.get_closed_loop_trajectory(
            MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
            concrete_controller,
            x0,
            nstep;
            stopping = reached,
        );
        println(control_trajectory)
        println()

        for k in 1:ST.length(control_trajectory)
            concrete_state = ST.get_state(control_trajectory, k)
            abstract_state = SY.get_abstract_state(abstract_system, concrete_state)
            println(abstract_state)
        end
        println()

        nstep = 300 # correspond to 30 sec
        function reached_abstract(x)
            if x ∈ abstract_problem.target_set
                return true
            else
                return false
            end
        end

        #######################################################
        ################# Abstract Trajectory #################
        #######################################################
        abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
        abstract_init_set = Dionysos.Symbolic.get_abstract_state(abstract_system, x0)
        #println(abstract_init_set)
        control_trajectory = get_abstract_closed_loop_trajectory(
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system")),
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller")),
            abstract_init_set,
            nstep;
            stopping = reached_abstract,
        );
        
        println(control_trajectory)

        #concrete_control_trajectory = get_concrete_trajectory(abstract_system, control_trajectory)
        #println(concrete_control_trajectory)

        println()
    end
    if(Second_part)
        println("Second Part")
        println()
        # Reload the result
        file = jldopen(filename_save, "r")
        reloaded_solver = file["my_abstraction_solver"]
        MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"), reloaded_solver)
        optimizer.handle_out_of_domain = 2

        #######################################################
        ################# Problem definition ##################
        #######################################################
        p0 = SVector{n_state,Float64}([-0.15352800685754736, 0.11944498327439435, 0.21311298746900986, 0.0, 0.0, 0.0])
        _I_ = UT.HyperRectangle(p0, p0)

        t_low = SVector{n_state,Float64}([-1.1*π/180.0, -1.1*π/180.0, -1.1*π/180.0, -0.75, -0.3, -0.3])
        t_high = SVector{n_state,Float64}([1.1*π/180.0, 1.1*π/180.0, 1.1*π/180.0, 0.3, 0.75, 0.75])
        _T_ = UT.HyperRectangle(t_low, t_high)
        
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
        abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
        concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
        controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
        uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))
        
        nstep = 300 # correspond to 30 sec
        function reached_abstract(x)
            if x ∈ abstract_problem.target_set
                return true
            else
                return false
            end
        end

        #######################################################
        ################# Abstract Trajectory #################
        #######################################################
        abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
        abstract_init_set = Dionysos.Symbolic.get_abstract_state(abstract_system, p0)
        #println(abstract_init_set)
        control_trajectory = get_abstract_closed_loop_trajectory(
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system")),
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller")),
            abstract_init_set,
            nstep;
            stopping = reached_abstract,
        );
        
        println(control_trajectory)

        concrete_control_trajectory = get_concrete_trajectory(abstract_system, control_trajectory)
        println(concrete_control_trajectory)

        println()

        
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
            MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
            concrete_controller,
            p0,
            nstep;
            stopping = reached,
            handle_out_of_domain = optimizer.handle_out_of_domain,
            state_space = state_space
        );

        println(concrete_control_trajectory)
        println()
        
    end
end