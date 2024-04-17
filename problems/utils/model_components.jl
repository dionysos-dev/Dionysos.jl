using ModelingToolkit, RigidBodyDynamics, DifferentialEquations

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

# function torque_ext(t, kt; exp = 1)
#     u_function = u_functions[exp]
#     voltage = u_function(t)
#     return voltage * kt
# end

# @component function VariableInput(; name, experiment = 1)
#     @named output = RealOutput()

#     u_function = u_functions[experiment]
#     equation = u_function(t)

#     eqs = [
#         output.u ~ equation
#     ]

#     compose(ODESystem(eqs, t, [], []; name = name), [output])
# end

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

@mtkmodel OLController begin
    @components begin
        input_voltage = RealInput()
        joint_out = Joint()
    end
    @parameters begin
        kt = 0.2843
        constant_input = false
    end
    @variables begin
        q(t) = 0., [description = "Absolute rotation angle", guess = 0.0]
        q̇(t) = 0., [description = "Absolute angular velocity", guess = 0.0]
    end
    @equations begin
        q ~ joint_out.q
        D(q) ~ q̇
        joint_out.τ ~ input_voltage.u * kt
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
        q(t) = 0., [description = "Absolute rotation angle", guess = 0.0]
        q̇(t) = 0., [description = "Absolute angular velocity", guess = 0.0]
        q̈(t) = 0., [description = "Absolute angular acceleration", guess = 0.0]
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
@connector function RealOutput(; name)
    @variables u(t)=0. [
        output = true,
        description = "Inner variable in RealOutput $name",
    ]
    ODESystem(Equation[], t, [u...], []; name = name)
end

@connector function RealInput(; name)
    @variables u(t)=0. [
        input = true,
        description = "Inner variable in RealInput $name"
    ]
    ODESystem(Equation[], t, [u...], []; name = name)
end

@mtkmodel PositionSensor begin
    @components begin
        joint_in = Joint()
        q = RealOutput()
    end
    @equations begin
        joint_in.q ~ q.u
        joint_in.τ ~ 0
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