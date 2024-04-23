using ModelingToolkit, RigidBodyDynamics, DifferentialEquations
using ModelingToolkitStandardLibrary.Blocks: RealOutput, RealInput, RealOutputArray, RealInputArray
### ---- Handling of the reference torque in Open-Loop ----

# Include using a relative path
# include("../utils/torque_references/references.jl")
# include("../utils/recorded_data/load_real_data.jl")

# # Using the module or its functions
# using .torque_references: waveform_value_at_time
# using .torque_references: current_for_time

# @register_symbolic (f::UFunction)(t)

# u_functions = []
# for i in eachindex(dfs)
#     push!(u_functions, UFunction(dfs[i], 2)) # Select between 3 methods
# end

# @component function VariableInput(; name, experiment = 1, nin = 1)
#     if nin == 1
#         @named output = RealOutput()
#         u_function = u_functions[experiment]
#         equation = u_function(t)
#     else
#         @named output = RealOutputArray(nout=nin)
#         equations = Array{typeof(u_functions[1](t))}(undef, nin)  # Initialize an array for the equations.
#         for i in 1:nin
#             u_function = u_functions[experiment[i]]  # Assume experiment is now an array of indices.
#             equations[i] = u_function(t)  # Apply each function to `t`.
#         end
#     end

#     eqs = [
#         output.u ~ equation
#     ]

#     compose(ODESystem(eqs, t, [], []; name = name), [output])
# end

# 22 April 2024 : MTK have not the possibility yet to create variable array length https://github.com/SciML/ModelingToolkit.jl/issues/2453
# So a quick and dirty bypass is proposed
@mtkmodel ConstantInput begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        U = 0.0, [description = "Constant output value"]
    end
    @equations begin
        output.u ~ U
    end
end

@mtkmodel ConstantInput2 begin
    @components begin
        output = RealOutputArray(nout=2)
    end
    @parameters begin
        U[1:2]::Float64 = zeros(2), [description = "Constant output value"]
    end
    @equations begin
        output.u[1] ~ U[1]
        output.u[2] ~ U[2]
    end
end


### ---- Compute the gravitational torque using RBD.jl ----

function compute_dynamics(mechanism, θs, ωs)
    actual_state = MechanismState(mechanism, θs, ωs)
    zero_velocity!(actual_state)
    v̇ = similar(velocity(actual_state))
    fill!(v̇, 0.0)
    τ_dynamics = -inverse_dynamics(actual_state, v̇)[1]
    return τ_dynamics
end

### ---- Components definition for robot modeling ----

@connector Joint begin
    q(t), [description = "Rotation angle of joint"]
    τ(t), [connect = Flow, description = "Cut torque in joint"]
end

# Same here, the @structural_parameters solution doesn't yet propose the possibility to define variable length parameters
# @mtkmodel OLController begin
#     @components begin
#         input_voltage = RealInput()
#         joint_out = Joint()
#     end
#     @parameters begin
#         kt = 0.2843
#         constant_input = false
#     end
#     @variables begin
#         q(t) = 0., [description = "Absolute rotation angle"]
#         q̇(t) = 0., [description = "Absolute angular velocity"]
#     end
#     @equations begin
#         q ~ joint_out.q
#         D(q) ~ q̇
#         joint_out.τ ~ input_voltage.u * kt
#     end
# end
@mtkmodel OLController begin
    @structural_parameters begin
        nin = 1
    end
    @components begin
        input_voltage = RealInput()
        joint_out = Joint()
        if nin == 2
            input_voltage = RealInputArray(nin=nin)
            joint_out = [Joint() for i in 1:2]
        end
    end
    @parameters begin
        kt = 0.2843
        constant_input = false
    end
    @variables begin
        q(t) = zeros(nin), [description = "Absolute rotation angle"]
        q̇(t) = zeros(nin), [description = "Absolute angular velocity"]
    end
    @equations begin
        if nin==1
            q[1] ~ joint_out.q
            D(q[1]) ~ q̇[1]
            joint_out.τ ~ input_voltage.u * kt
        else
            for i in 1:nin
                q[2] ~ joint_out[i].q
                D(q[2]) ~ q̇[2]
                joint_out[i].τ ~ input_voltage[i].u * kt
            end
        end
    end
end

@mtkmodel OLController2 begin
    @components begin
        input_voltage = RealInputArray(nin=2)
        # joint_out = [Joint() for i in 1:2] -> doesn't work...
        joint_out1 = Joint()
        joint_out2 = Joint()
    end
    @parameters begin
        kt = 0.2843
        constant_input = false
    end
    @variables begin
        q(t)[1:2] = zeros(2), [description = "Absolute rotation angle"]
        q̇(t)[1:2] = zeros(2), [description = "Absolute angular velocity"]
    end
    @equations begin
        D(q) ~ q̇
        q[1] ~ joint_out1.q
        joint_out1.τ ~ input_voltage.u[1] * kt
        q[2] ~ joint_out2.q
        joint_out2.τ ~ input_voltage.u[2] * kt
    end
end

@mtkmodel GravitationalTorque begin
    @structural_parameters begin
        hip = false
        single_pendulum = true
    end
    @components begin
        joint_in = Joint()
        joint_out = Joint()
        if !single_pendulum 
            other_joint = Joint()
        end
    end
    @variables begin
        q(t), [description = "Absolute rotation angle", guess = 0.0]
        q̇(t), [description = "Absolute angular velocity", guess = 0.0]
    end
    begin
        if single_pendulum
            τ_grav = compute_dynamics(mechanism, [q], [q])
        elseif hip
            τ_grav = compute_dynamics(mechanism, [q, other_joint.q], [q, other_joint.q])
        else # knee of double_pendulum
            τ_grav = compute_dynamics(mechanism, [other_joint.q, q], [other_joint.q, q])
        end
    end
    @equations begin
        q ~ joint_in.q
        q ~ joint_out.q
        D(q) ~ q̇
        joint_out.τ ~ -joint_in.τ + τ_grav
    end
end

@mtkmodel FrictionTorque begin
    @structural_parameters begin
        stribeck = false
    end
    @components begin
        joint_in = Joint()
        joint_out = Joint()
    end
    @parameters begin
        τ_c = 0.128, [description = "Coulomb torque"]
        Kv = 1.2465, [description = "Viscous coefficient"]
        τ_s = 0.4, [description = "Static friction torque"]
        q̇_s = 0.1, [description = "Stribeck coefficient"]
    end
    @variables begin
        q(t), [description = "Absolute rotation angle", guess = 0.0]
        q̇(t), [description = "Absolute angular velocity", guess = 0.0]
    end
    begin
        if stribeck
            str_scale = (τ_s - τ_c)
            β = exp(-abs(q̇/q̇_s))
            τ_f = (τ_c - str_scale * β) * tanh(q̇ / 1e-4) + Kv * q̇
        else
            τ_f = tanh(q̇ / 1e-4) * τ_c + Kv * q̇
        end
    end
    @equations begin
        q ~ joint_in.q
        q ~ joint_out.q
        D(q) ~ q̇
        joint_out.τ ~ -joint_in.τ - τ_f
    end
end

@mtkmodel TotalInertia begin
    @structural_parameters begin
        hip = true
        single_pendulum = true
    end
    @parameters begin
        J_motor = 0.07305, [description = "Motor's Moment of inertia"]
    end
    @components begin
        joint_in = Joint()
        joint_out = Joint()
        if hip 
            joint_knee = Joint()
        end
    end
    # begin
    #     @symcheck J_motor > 0 || throw(ArgumentError("Expected `J` to be positive"))
    # end
    @variables begin
        q(t) = 0., [description = "Absolute rotation angle"]
        q̇(t) = 0., [description = "Absolute angular velocity"]
        q̈(t) = 0., [description = "Absolute angular acceleration"]
    end
    begin
        if single_pendulum
            J_external = mass_matrix(MechanismState(mechanism, [q], [q̇]))[1]
        elseif hip
            J_external = mass_matrix(MechanismState(mechanism, [q, joint_knee.q], [q, joint_knee.q]))[1]
        else 
            J_external = mass_matrix(MechanismState(mechanism, [0., q], [0., q]))[4]
        end
    end
    @equations begin
        q ~ joint_in.q
        q ~ joint_out.q
        D(q) ~ q̇
        D(q̇) ~ q̈
        q̈ ~ -joint_in.τ / (J_motor + J_external)
    end
end

### ---- Components and connectors definition needed for linearization (source : ModelingToolkitStandardLibrary.jl) ----
@mtkmodel PositionSensor begin
    @structural_parameters begin
        nin = 1
    end
    @components begin
        joint_in = Joint()
        q = RealOutput()
        if nin > 1
            joint_in = Joint()
            q = RealOutputArray(nout = nin)
        end
    end
    @equations begin
        if nin == 1
            joint_in.q ~ q.u
            joint_in.τ ~ 0
        else 
            joint_in.q ~ q.u
            joint_in.τ ~ 0
        end
    end
end

@mtkmodel PositionSensor2 begin
    @components begin
        joint_in1 = Joint()
        joint_in2 = Joint()
        q = RealOutputArray(nout = 2)
    end
    @equations begin
        joint_in1.q ~ q.u[1]
        joint_in1.τ ~ 0
        joint_in2.q ~ q.u[2]
        joint_in2.τ ~ 0
    end
end

@mtkmodel Feedback1D begin
    @components begin
        input = RealInput()
        output = RealOutput()
    end
    @equations begin
        output.u ~ input.u
    end
end

@mtkmodel Feedback2D begin
    @components begin
        input = RealInputArray(nin = 2)
        output = RealOutputArray(nout = 2)
    end
    @equations begin
        output.u ~ input.u
    end
end