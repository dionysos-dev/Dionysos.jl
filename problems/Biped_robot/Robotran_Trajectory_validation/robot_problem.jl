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

using Libdl
# Load library
lib = Libdl.dlopen(joinpath(@__DIR__, "../../../../philippides_J2C/workR/build/libProject_user.so"))
philippides_func = Libdl.dlsym(lib, :philippides)
get_res_func    = Libdl.dlsym(lib, :get_philippides_results)
function call_philippides(x::Vector{Float64})
    ccall(philippides_func, Cvoid, (Ptr{Float64},), x)
end
function get_results()
    res = Vector{Float64}(undef, 8)
    ccall(get_res_func, Cvoid, (Ptr{Float64},), res)
    return res
end

# include the tools for the simulator from src

function system(;
    tstep = 5e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    ## Robots Parameters ##
    Lthigh = 0.20125
    Lleg = 0.172
    Hip_offset = 0.04025
    Foot_height = 0.009
    Init_offset = 0.0 # Was only valid when using Juliarobotics

    function fill_state!(x)
        # Create q
        q = vcat(zeros(2), x[1:3], zeros(3))
        q̇ = vcat(zeros(2), x[4:6], zeros(3))

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
        q_ref = SVector{1,Float64}(0.0)
        results = []
        cd(joinpath(@__DIR__,"../../../../philippides_J2C/workR/build")) do
            x = [q...,q̇...,u...,q_ref...]
            call_philippides(x)
            res = get_results()
            push!(results,res...)
        end        

        x_next = SVector{6}(results[1:3]...,results[5:7]...)
        return x_next
    end
    # Define state space (bounds should be set according to your robot's joint limits)
    # Note : We need to add the discretisation step at each of the borns if the ones we chose are supposed to be centroids
    disc_steps = [fill(π/180, 3)..., fill(0.075, 3)...]
    state_lower_bounds = [-12*π/180, 0, 0, -0.6, -0.3, -0.6] .- disc_steps
    state_upper_bounds = [0, 12*π/180, 14*π/180, 0.3, 0.6, 0.6] .+ disc_steps

    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    # Define input space (bounds should be set according to actuator limits)
    # Note : We con't need for the inputs to add something if the born are centroids
    input_lower_bounds = [-3, -3, -3]   # Example: torque or force limits
    input_upper_bounds = [3, 2, 3]

    input_space = UT.HyperRectangle(input_lower_bounds, input_upper_bounds)

    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        3 + 3, # state space : the 3 actuators in position and speed (right knee not included)
        3,     # input space : the volatge on 3 actuators (right knee not included)
        state_space,
        input_space,
    )

    return sys
end

function problem(;
    tstep = 1e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf)
    return PB.EmptyProblem(sys, nothing)
end

end
