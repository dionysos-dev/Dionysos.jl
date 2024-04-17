module double_pendulum_system

### Please note : the sign(velocity) function is a tanh, because continuous/discrete events are not (not yet ?) working well with components (I didn't achieve to make it working)

using ModelingToolkit, Plots, DataFrames, RigidBodyDynamics, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D

include("../model_components.jl")

function system(mechanism_; pars = [0.128, 1.2465, 0.2843, 0.07305, 0.128, 0], stribeck=false, experiment = 1, constant_input=false)
    global mechanism = mechanism_
    if !stribeck
        push!(pars, 0.)
    end
    if constant_input
        @named v_input = ConstantInput(U = experiment[1])
        @named v_input2 = ConstantInput(U = experiment[2])
    else
        @named v_input = VariableInput(experiment = experiment[1])
        @named v_input2 = VariableInput(experiment = experiment[2])
    end

    # Upper link
    @named ctrl = OLController(kt = pars[3])
    @named J_total = TotalInertia(J_motor = pars[4], single_pendulum = false)
    @named τ_f = FrictionTorque(stribeck = stribeck, τ_c = pars[1], Kv = pars[2], τ_s = pars[5], q̇_s = pars[6])
    @named τ_grav = GravitationalTorque(hip = true, single_pendulum = false)
    @named position_sensor = PositionSensor()
    @named output_reader = Feedback1D()

    # Lower link
    @named ctrl2 = OLController(kt = pars[3])
    @named J_total2 = TotalInertia(J_motor = pars[4], hip = false, single_pendulum = false)
    @named τ_f2 = FrictionTorque(stribeck = stribeck, τ_c = pars[1], Kv = pars[2], τ_s = pars[5], q̇_s = pars[6])
    @named τ_grav2 = GravitationalTorque(single_pendulum = false)
    @named position_sensor2 = PositionSensor()
    @named output_reader2 = Feedback1D()

    connections = [ connect(v_input.output, :i, ctrl.input_voltage)
                    connect(ctrl.joint_out, τ_grav.joint_in)
                    connect(J_total2.joint_out, τ_grav.other_joint)
                    connect(τ_grav.joint_out, τ_f.joint_in)
                    connect(τ_f.joint_out,  J_total.joint_in)
                    connect(J_total2.joint_out, J_total.joint_knee)
                    connect(J_total.joint_out, position_sensor.joint_in)
                    connect(position_sensor.q, :o, output_reader.input)

                    connect(v_input2.output, :i2, ctrl2.input_voltage)
                    connect(ctrl2.joint_out, τ_grav2.joint_in)
                    connect(J_total.joint_out, τ_grav2.other_joint)
                    connect(τ_grav2.joint_out, τ_f2.joint_in)
                    connect(τ_f2.joint_out,  J_total2.joint_in)
                    connect(J_total2.joint_out, position_sensor2.joint_in)
                    connect(position_sensor2.q, :o2, output_reader2.input)
                    ]

    @named model = ODESystem(connections, t,
                        systems = [
                            v_input,
                            ctrl,
                            J_total,
                            τ_f,
                            τ_grav,
                            position_sensor,
                            output_reader,
                            v_input2,
                            ctrl2,
                            J_total2,
                            τ_f2,
                            τ_grav2,
                            position_sensor2,
                            output_reader2
                        ])
    sys = structural_simplify(model)
    return sys, model
end

function problem(sys; u0 = [0., 0., 0., 0.], tspan = (0.0, 10.0), constant_input = false, experiment = 1)
    if !constant_input
        tspan = (dfs[experiment[1]].timestamp[1], dfs[experiment[1]].timestamp[end])
    end
    prob = ODEProblem(sys, [sys.J_total.q => u0[1], sys.J_total.q̇ => u0[2], sys.J_total2.q => u0[3], sys.J_total2.q̇ => u0[4]], tspan, [])
    return prob
end

function get_next_state(sys, uk, input, tspan)
    prob = ODEProblem(sys, [sys.J_total.q => uk[1], sys.J_total.q̇ => uk[2]], (0., tspan), [sys.v_input.U => input])
    sol = solve(prob, Rosenbrock23(), dtmax=0.01, reltol=1e-3, abstol=1e-3)
    next_u = (sol.u[end][1], sol.u[end][2])
    return next_u
end

function linear_system(sys, model, operation_point)
    op = Dict(sys.J_total.q => operation_point[1], sys.J_total.q̇ => operation_point[2])
    matrices, simplified_sys = linearize(model, :i, :o, op=op)
    return matrices
end

end