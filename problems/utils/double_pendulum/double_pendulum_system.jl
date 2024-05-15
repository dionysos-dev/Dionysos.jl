module double_pendulum_system

### Please note : the sign(velocity) function is a tanh, because continuous/discrete events are not (not yet ?) working well with components (I didn't achieve to make it working)

using ModelingToolkit, RigidBodyDynamics, DifferentialEquations, ModelingToolkitStandardLibrary
using ModelingToolkit: t_nounits as t, D_nounits as D

include("../model_components.jl")

function system(mechanism_; pars = [0.06510461345450586, 1.5879662676966781, 0.39454422423683916, 0., 0.06510461345450586, 0.], stribeck=false, experiment = nothing, constant_input=[0., 0.])
    global mechanism = mechanism_
    if !stribeck
        push!(pars, 0.)
    end
    if experiment !== nothing
        @named v_input = VariableInput2(experiment = experiment, nin = 2)
    else
        @named v_input = ConstantInput2(U = constant_input)
    end

    @named ctrl = OLController2(kt = pars[3])
    @named position_sensor = PositionSensor2()
    @named output_reader = Feedback2D()

    # Upper link
    @named J_total = TotalInertia(J_motor = pars[4], single_pendulum = false)
    @named τ_f = FrictionTorque(stribeck = stribeck, τ_c = pars[1], Kv = pars[2], τ_s = pars[5], q̇_s = pars[6])
    @named τ_grav = GravitationalTorque(hip = true, single_pendulum = false)

    # Lower link
    @named J_total2 = TotalInertia(J_motor = pars[4], hip = false, single_pendulum = false)
    @named τ_f2 = FrictionTorque(stribeck = stribeck, τ_c = pars[1], Kv = pars[2], τ_s = pars[5], q̇_s = pars[6])
    @named τ_grav2 = GravitationalTorque(single_pendulum = false)

    connections = [ connect(v_input.output, :i, ctrl.input_voltage)

                    connect(ctrl.joint_out1, τ_grav.joint_in)
                    connect(J_total2.joint_out, τ_grav.other_joint)
                    connect(τ_grav.joint_out, τ_f.joint_in)
                    connect(τ_f.joint_out,  J_total.joint_in)
                    connect(J_total2.joint_out, J_total.joint_knee)
                    connect(J_total.joint_out, position_sensor.joint_in1)

                    connect(ctrl.joint_out2, τ_grav2.joint_in)
                    connect(J_total.joint_out, τ_grav2.other_joint)
                    connect(τ_grav2.joint_out, τ_f2.joint_in)
                    connect(τ_f2.joint_out,  J_total2.joint_in)
                    connect(J_total2.joint_out, position_sensor.joint_in2)

                    connect(position_sensor.q, :o, output_reader.input)
                    ]

    @named model = ODESystem(connections, t,
                        systems = [
                            v_input,
                            ctrl,
                            J_total,
                            τ_f,
                            τ_grav,
                            J_total2,
                            τ_f2,
                            τ_grav2,
                            position_sensor,
                            output_reader,
                        ])
    sys = structural_simplify(model)
    return sys, model
end

function problem(sys; u0 = [0., 0., 0., 0.], tspan = (0.0, 10.0), experiment = nothing, constant_input=nothing)
    if experiment !== nothing
        tspan = (dfs[experiment[1]].timestamp[1], dfs[experiment[1]].timestamp[end])
        prob = ODEProblem(sys, [sys.J_total.q => u0[1], sys.J_total.q̇ => u0[2], sys.J_total2.q => u0[3], sys.J_total2.q̇ => u0[4]], tspan, [])
    else
        prob = ODEProblem(sys, [sys.J_total.q => u0[1], sys.J_total.q̇ => u0[2], sys.J_total2.q => u0[3], sys.J_total2.q̇ => u0[4]], tspan, [sys.v_input.U[1] => constant_input[1], sys.v_input.U[2] => constant_input[2]])
    end
    return prob
end

function get_next_state(sys, prob, uk, input, tspan)
    newprob = remake(prob, tspan = (0, tspan); u0 = [sys.J_total.q => uk[1], sys.J_total.q̇ => uk[2], sys.J_total2.q => uk[3], sys.J_total2.q̇ => uk[4]], p = [sys.v_input.U[1] => input[1], sys.v_input.U[2] => input[2]])
    sol = solve(newprob, Tsit5(), dtmax=0.01, reltol=1e-3, abstol=1e-3)
    next_u = (sol.u[end][1], sol.u[end][2], sol.u[end][3], sol.u[end][4])
    return next_u
end

function linear_system(sys, model, operation_point)
    op = Dict(sys.J_total.q => operation_point[1], sys.J_total.q̇ => operation_point[2], sys.J_total2.q => operation_point[3], sys.J_total2.q̇ => operation_point[4])
    matrices, simplified_sys = linearize(model, :i, :o, op=op)
    return matrices
end

end