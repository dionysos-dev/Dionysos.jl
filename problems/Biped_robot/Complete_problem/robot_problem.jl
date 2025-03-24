module RobotProblem

using MathematicalSystems
using StaticArrays
using RigidBodyDynamics

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic

# include the tools for the simulator from src
include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

function get_visualization_tool(;
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    # Construct the robot in the simulation engine 
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true, # TODO ->true
    )
    vis = RS_tools.set_visulalizer(; mechanism = rs.mechanism, fileName = robot_urdf)
    return rs, vis
end

function system(;
    tstep = 5e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )

    mechanism = rs.mechanism
    state = MechanismState(mechanism)
    n_pos = num_positions(state)
    n_vel = num_velocities(state)
    Δt_simu = 1e-4       # Simulation step 
    Δt_dionysos = tstep      # Dinoysos time discretisation, nominal 50Hz (control freq of the material robot)

    println("n_pos: ", n_pos)
    println("n_vel: ", n_vel)

    ## MOTOR Parameters ##
    HGR = 353.5           # Hip gear-ratio
    KGR = 212.6           # Knee gear-ratio
    ktp = 0.395 / HGR      # Torque constant with respect to the voltage [Nm/V] 
    Kvp = 1.589 / (HGR * HGR)      # Viscous friction constant [Nm*s/rad] (linked to motor speed)
    τc_u = 0.065 / HGR           # Dry friction torque [Nm]
    τ_m = [0.0, 0.0, 0.0, 0.0]
    # Discrete time using Rigibodydynamics simulator -> returns (X[i], U[i]) -> X[i+1]
    function voltage_controller!(u::SVector)
        ddl = 2
        function controller!(τ, t, state)
            τ .= 0
            current_̇q = velocity(state)[(end - 3 - ddl):(end - ddl)]
            ω = current_̇q .* [HGR, HGR, KGR, KGR]

            τ_0 = u .* [HGR, HGR, KGR, KGR] .* ktp .- ω .* [HGR, HGR, KGR, KGR] .* Kvp
            τ_m .= τ_0 .- sign.(ω) .* [HGR, HGR, KGR, KGR] .* τc_u
            τ[(end - 3 - ddl):(end - ddl)] .= τ_m
            return nothing
        end
    end

    ## Robots Parameters ##
    Lthigh = 0.20125
    Lleg = 0.172
    Hip_offset = 0.04025
    Foot_height = 0.009
    Init_offset = -0.0006559432

    function fill_state!(x)
        # Create q
        q = vcat(zeros(2), x[1:4], zeros(2))
        q̇ = vcat(zeros(2), x[5:8], zeros(2))

        # Compute the heights of the two legs (double pendulums)
        zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
        zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])

        # FILL THE POSITIONS

        # Write the maximum height to q[2]
        # (adding the distance from the hip joint to hip body and the height of the foot)
        # The most extended leg is in contact with the ground
        
        # The height of the boom is set at 0 ! # Note: there is a slight error in the URDF and the robot is flying => Init_offset
        q[2] = max(zl, zr) - Lthigh - Lleg + Init_offset

        # Set additional constraints
        # The x position is set to 0
        # The feet are kept // to the ground 

        q[7] = -(q[3] + q[5])
        q[8] = -(q[4] + q[6])

        # FILL THE SPEEDS
        # identify the contact leg
        i1 = 0
        i2 = 0
        if (zl > zr)
            i1, i2 = 3, 5
        else
            i1, i2 = 4, 6
        end
        # speed equations of the double pendulum
        x = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
        ẋ = Lthigh * q̇[i1] * cos(q[i1]) + Lleg * (q̇[i1] + q̇[i2]) * cos(q[i1] + q[i2])
        ż = -(Lthigh * q̇[i1] * sin(q[i1]) + Lleg * (q̇[i1] + q̇[i2]) * sin(q[i1] + q[i2]))
        q[1] = x
        q̇[1] = ẋ
        q̇[2] = ż

        # adjust the angular speed of the feet to remain mostly horizontal
        q̇[7] = -(q̇[3] + q̇[5])
        q̇[8] = -(q̇[4] + q̇[6])
        return q, q̇
    end

    function vectorFieldBipedRobot(x, u)
        # Variables: [x z LH RH LK RK LA RA]
        # NB: to move the knee forward, a negative angle is needed!

        # First step: fill state: from the n state variables -> 8 positions and 8 speeds
        q, q̇ = fill_state!(x)
        
        # Second step: set the mechanism in that configuration
        set_configuration!(state, q)
        set_velocity!(state, q̇)

        # Third step: get next state
        controller! = voltage_controller!(u)
        ts, qs, vs =
            RigidBodyDynamics.simulate(state, Δt_dionysos, controller!; Δt = Δt_simu)
        x_next = SVector{length(x)}(qs[end][3:6]..., vs[end][3:6]...) 
        # Note: qs and vs are vectors of speed and position for every step of the simulation (i.e. every Δt = 1e-4)
        # Only the final states are useful in our case

        return x_next
    end

    # Define state space (bounds should be set according to your robot's joint limits)
    state_lower_bounds = [-0.5, -0.5, -0.2, -0.2, -0.8, -0.8, -0.8, -0.8]  # W350 -> max 30rpm, W210 -> max 50rpm
    state_upper_bounds = [0.5 ,  0.5,  0.8,  0.8,  0.8,  0.8,  0.8,  0.8]

    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    # Define input space (bounds should be set according to actuator limits)
    input_lower_bounds = [-3, -3, -3, -3]   # Example: torque or force limits
    input_upper_bounds = [3, 3, 3, 3]

    input_space = UT.HyperRectangle(input_lower_bounds, input_upper_bounds)

    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        4 + 4, # state space : the 4 actuators in position and speed
        4,     # input space : the volatge on 4 actuators
        state_space,
        input_space,
    )

    return sys
end

function problem(;
    tstep = 2e-2,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf)
    return PB.EmptyProblem(sys, nothing)
end

function deep_equal(a, b)
    # If both are the same reference, they are equal
    if a === b
        return true
    end

    # Handle primitive types and strings directly
    if isbits(a) || isa(a, AbstractString)
        return a == b
    end

    # Handle arrays element-wise
    if isa(a, AbstractArray)
        return length(a) == length(b) && all(deep_equal(a[i], b[i]) for i in eachindex(a, b))
    end

    # Handle dictionaries by comparing keys and values
    if isa(a, AbstractDict) && isa(b, AbstractDict)
        return length(a) == length(b) && all(k -> deep_equal(a[k], b[k]), keys(a))
    end

    # Handle Set elements by checking equality, disregarding order
    if isa(a, Set) && isa(b, Set)
        return length(a) == length(b) && all(x -> x in b, a) && all(x -> x in a, b)
    end

    # Handle Any-type fields recursively (if both are Any)
    if isa(a, Any) || isa(b, Any)
        return deep_equal(a, b)
    end

    # Handle Grid-like structures (e.g., DomainList, GridFree)
    if isa(a, GridDomainType) && isa(b, GridDomainType)
        return deep_equal(a.grid, b.grid) && deep_equal(a.elems, b.elems)
    end

    # Recursively compare fields for composite types (structs, tuples, named tuples)
    for field in fieldnames(typeof(a))
        val_a = getfield(a, field)
        val_b = getfield(b, field)
        if !deep_equal(val_a, val_b)
            println("Difference found in field: ", field)
            return false
        end
    end

    return true
end

end
