module single_pendulum_system

### Please note : the sign(velocity) function is a tanh, because continuous/discrete events are not (not yet ?) working well with components (I didn't achieve to make it working)

using ModelingToolkit, ModelingToolkitStandardLibrary, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D

include(joinpath("..", "model_components.jl"))

function system(mechanism_; pars = [0.06510461345450586, 1.5879662676966781, 0.39454422423683916, 0., 0.06510461345450586, 0.], stribeck=false, experiment = nothing, constant_input = 0.)
    global mechanism = mechanism_
    if experiment !== nothing
        @named v_input = VariableInput(experiment = experiment)
    else
        @named v_input = ConstantInput(U = constant_input)
    end
    @named ctrl = OLController(kt = pars[3])
    @named J_total = TotalInertia(J_motor = pars[4])
    @named τ_f = FrictionTorque(stribeck = stribeck, τ_c = pars[1], Kv = pars[2], τ_s = pars[5], q̇_s = pars[6])
    @named τ_grav = GravitationalTorque()
    @named position_sensor = PositionSensor()
    @named output_reader = Feedback1D()

    connections = [ connect(v_input.output, :i, ctrl.input_voltage)
                    connect(ctrl.joint_out, τ_grav.joint_in)
                    connect(τ_grav.joint_out, τ_f.joint_in)
                    connect(τ_f.joint_out,  J_total.joint_in)
                    connect(J_total.joint_out, position_sensor.joint_in)
                    connect(position_sensor.q, :o, output_reader.input)
                    ]

    @named model = ODESystem(connections, t,
                        systems = [
                            v_input,
                            ctrl,
                            J_total,
                            τ_f,
                            τ_grav,
                            position_sensor,
                            output_reader
                        ])
    sys = structural_simplify(model)
    return sys, model
end

function problem(sys; u0 = [0., 0.], tspan = (0.0, 10.0), experiment = nothing, constant_input=nothing)
    if experiment !== nothing
        tspan = (dfs[experiment].timestamp[1], dfs[experiment].timestamp[end])
        prob = ODEProblem(sys, [sys.J_total.q => u0[1], sys.J_total.q̇ => u0[2]], tspan, [])
    else 
        prob = ODEProblem(sys, [sys.J_total.q => u0[1], sys.J_total.q̇ => u0[2]], tspan, [sys.v_input.U => constant_input])
    end
    return prob
end

function get_next_state(sys, prob, uk, input, tspan)
    newprob = remake(prob, tspan = (0, tspan); u0 = [sys.J_total.q => uk[1], sys.J_total.q̇ => uk[2]], p = [sys.v_input.U => input])
    sol = solve(newprob, Tsit5(), dtmax=0.01, reltol=1e-3, abstol=1e-3)
    next_u = (sol.u[end][1], sol.u[end][2])
    return next_u
end

function linear_system(sys, model, operation_point)
    op = Dict(sys.J_total.q => operation_point[1], sys.J_total.q̇ => operation_point[2])
    matrices, simplified_sys = linearize(model, :i, :o, op=op)
    return matrices
end

end