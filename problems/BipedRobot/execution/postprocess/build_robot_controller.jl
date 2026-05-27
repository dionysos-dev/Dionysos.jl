using Printf
using Dionysos
using MathematicalSystems
using StaticArrays
using LinearAlgebra
using JuMP
using MathOptInterface
using JLD2

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const SY = DI.Symbolic
const MP = DI.Mapping
const ST = DI.System

const TSTEP = 0.1

const ROBOT_PROBLEM_FILE = joinpath(@__DIR__, "..", "..", "6D_model", "robot_problem.jl")
const ROBOT_URDF_FILE =
    joinpath(@__DIR__, "..", "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")
const ABSTRACTION_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_abstraction_1.jld2")
const CONTROLLER_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_controller_2.jld2")

include(ROBOT_PROBLEM_FILE)
using .RobotProblem

function load_optimizer(filename::AbstractString)
    return JLD2.jldopen(filename, "r") do file
        return file["optimizer"]
    end
end

function save_controller(
    filename::AbstractString;
    controller,
    concrete_system,
    control_problem,
)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        file["controller"] = controller
        file["concrete_system"] = concrete_system
        file["control_problem"] = control_problem
        return file["tstep"] = TSTEP
    end
    return filename
end

@info "Loading abstraction" ABSTRACTION_FILE
optimizer = load_optimizer(ABSTRACTION_FILE)

@info "Building concrete robot problem"
concrete_problem = RobotProblem.problem(; robot_urdf = ROBOT_URDF_FILE, tstep = TSTEP)
concrete_system = concrete_problem.system
n_state = MathematicalSystems.statedim(concrete_system)

x0 = SVector{n_state, Float64}(zeros(n_state))

target_low =
    SVector{n_state, Float64}([-12π / 180, 7π / 180, 8π / 180, -0.75, -0.30, -0.30])

target_high = SVector{n_state, Float64}([-8π / 180, 9π / 180, 12π / 180, 0.30, 0.75, 0.75])

I = UT.HyperRectangle(x0, x0)
T = UT.HyperRectangle(target_low, target_high)

control_problem = DI.Problem.OptimalControlProblem(
    concrete_system,
    I,
    T,
    nothing,
    nothing,
    DI.Problem.Infinity(),
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), control_problem)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("out_of_domain_handler"),
    SY.ProjectToNearestCellHandler(; warn = true),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

@info "Solving controller"
t_solve = @elapsed MOI.optimize!(optimizer)
@printf("Control solve time: %.3f s\n", t_solve)

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

save_controller(
    CONTROLLER_FILE;
    controller = controller,
    concrete_system = concrete_system,
    control_problem = control_problem,
)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
Xmap = SY.get_state_mapping(abstract_system)
println(
    "Number of controllable states = ",
    length(collect(MP.enum_states(controllable_set, Xmap))),
)

@info "Saved controller" CONTROLLER_FILE
