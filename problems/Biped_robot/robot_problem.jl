module RobotProblem

using MathematicalSystems
using LinearAlgebra, StaticArrays
using RigidBodyDynamics

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic

# include BipedRobot
include(joinpath(@__DIR__, "..", "..", "BipedRobot", "src", "BipedRobot.jl"))
using .BipedRobot

# include the tools for the simulator from src
include(joinpath(@__DIR__, "src", "RS_tools.jl"))
import .RS_tools

function get_visualization_tool(; robot_urdf = joinpath(@__DIR__, "deps/ZMP_2DBipedRobot_nodamping.urdf"))
    # Construct the robot in the simulation engine 
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    );

    vis = RS_tools.set_visulalizer(; mechanism = rs.mechanism, fileName = robot_urdf)
    return rs, vis
end


function system(; tstep = 2e-2, robot_urdf = joinpath(@__DIR__, "deps/ZMP_2DBipedRobot_nodamping.urdf"))
    robot = BipedRobot.mechanism(;
    symbolic = false,
    add_contact_points = true,
    add_flat_ground = true,)
    state = MechanismState(robot)
    n_pos = num_positions(state)
    n_vel = num_velocities(state)
    n_add = num_additional_states(state)
    Δt_simu     = 1e-4       # Simulation step 
    Δt_dionysos = tstep      # Dinoysos time discretisation, nominal 50Hz (control freq of the material robot)

    println("n_pos: ", n_pos)
    println("n_vel: ", n_vel)
    println("n_add: ", n_add)
    
    # Discrete time using Rigibodydynamics simulator -> returns (X[i], U[i]) -> X[i+1]
    function get_constant_controller!(
        u:: SVector,
    )
        sim_index = 0
        ddl = 2
        function controller!(τ, t, state)
            if(sim_index == 0)
                τ .= 0
                τ[(end - 3 - ddl):(end - ddl)] .= u
                sim_index = 1
            end
        end
    end

    function vectorFieldBipedRobot(x, u)
        q = x[1:n_pos]
        v = x[(n_pos + 1):(n_pos + n_vel)]
        s = x[(n_pos + n_vel + 1):end] # Utile ??
        println(length(x))
        set_configuration!(state, q)
        set_velocity!(state, v)
        set_additional_state!(state, s) # Utile ??
    
        controller! = get_constant_controller!(u)
        ts, qs, vs  = RigidBodyDynamics.simulate(state, Δt_dionysos, controller!; Δt = Δt_simu);
        
        x_next = SVector{2 * length(qs[end])}(qs[end]..., vs[end]...)
    
        return x_next
    end


    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        n_pos + n_vel + n_add,
        n_vel,
        nothing,
        nothing
    )

    return sys
end

function problem(; tstep = 2e-2, robot_urdf = joinpath(@__DIR__, "deps/ZMP_2DBipedRobot_nodamping.urdf"))
    println("LOAD Problem")
    sys = system(; tstep = tstep, robot_urdf = robot_urdf)
    return PB.EmptyProblem(sys, nothing)
end

end
