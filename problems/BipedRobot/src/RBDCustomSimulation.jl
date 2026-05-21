module RBDCustomSimulation

using RigidBodyDynamics

const RBD = RigidBodyDynamics
const OI = RigidBodyDynamics.OdeIntegrators

export simulate_final_state!,
    simulate_final_state_history!,
    simulate_final_state_custom!,
    CachedFinalStateSimulator,
    simulate_final_state_cached!

struct NoStorage{T} <: OI.OdeResultsSink end

Base.eltype(::NoStorage{T}) where {T} = T

OI.initialize(::NoStorage, t, state) = nothing
OI.process(::NoStorage, t, state) = nothing

function integrate_fixed!(integrator, final_time, Δt)
    T = eltype(integrator)
    t = zero(T)

    nsteps_float = final_time / Δt
    nsteps = round(Int, nsteps_float)

    if !isapprox(nsteps * Δt, final_time; rtol = 1e-10, atol = 1e-12)
        error(
            "final_time must be an integer multiple of Δt. Got final_time=$final_time, Δt=$Δt",
        )
    end

    OI.initialize(integrator.sink, t, integrator.state)

    for _ in 1:nsteps
        OI.step(integrator, t, Δt)
        t += Δt
        OI.process(integrator.sink, t, integrator.state)
    end

    return nothing
end

function simulate_final_state!(
    state::MechanismState,
    final_time,
    control! = RBD.zero_torque!;
    Δt,
    backend::Symbol = :history,
    kwargs...,
)
    if backend == :history
        return simulate_final_state_history!(
            state,
            final_time,
            control!;
            Δt = Δt,
            kwargs...,
        )
    elseif backend == :custom
        return simulate_final_state_custom!(state, final_time, control!; Δt = Δt, kwargs...)
    else
        error("Unknown simulation backend: $backend. Expected :history or :custom.")
    end
end

function simulate_final_state_history!(
    state::MechanismState,
    final_time,
    control! = RBD.zero_torque!;
    Δt,
    kwargs...,
)
    _, qs, vs = RBD.simulate(state, final_time, control!; Δt = Δt, kwargs...)

    return copy(qs[end]), copy(vs[end])
end

function simulate_final_state_custom!(
    state0::MechanismState{X},
    final_time,
    control! = RBD.zero_torque!;
    Δt = 1e-4,
    stabilization_gains = RBD.default_constraint_stabilization_gains(X),
    check_finite::Bool = true,
) where {X}
    sim = CachedFinalStateSimulator(
        state0,
        control!;
        stabilization_gains = stabilization_gains,
    )

    return simulate_final_state_cached!(sim, final_time, Δt; check_finite = check_finite)
end

mutable struct CachedFinalStateSimulator{C, I}
    controller!::C
    integrator::I
end

function CachedFinalStateSimulator(
    state0::MechanismState{X},
    controller!;
    stabilization_gains = RBD.default_constraint_stabilization_gains(X),
) where {X}
    T = RBD.cache_eltype(state0)

    result = RBD.DynamicsResult{T}(state0.mechanism)
    control_torques = similar(RBD.velocity(state0))

    closed_loop_dynamics! =
        let result = result,
            control_torques = control_torques,
            controller! = controller!,
            stabilization_gains = stabilization_gains

            function (v̇::AbstractArray, ṡ::AbstractArray, t, state)
                controller!(control_torques, t, state)

                RBD.dynamics!(
                    result,
                    state,
                    control_torques;
                    stabilization_gains = stabilization_gains,
                )

                copyto!(v̇, result.v̇)
                copyto!(ṡ, result.ṡ)

                return nothing
            end
        end

    tableau = OI.runge_kutta_4(T)
    storage = NoStorage{T}()

    integrator = OI.MuntheKaasIntegrator(state0, closed_loop_dynamics!, tableau, storage)

    return CachedFinalStateSimulator(controller!, integrator)
end

function simulate_final_state_cached!(
    sim::CachedFinalStateSimulator,
    final_time,
    Δt;
    check_finite::Bool = true,
)
    integrate_fixed!(sim.integrator, final_time, Δt)

    q_end = copy(RBD.configuration(sim.integrator.state))
    v_end = copy(RBD.velocity(sim.integrator.state))

    if check_finite
        if any(!isfinite, q_end) || any(!isfinite, v_end)
            @show final_time Δt q_end v_end
            error("Cached custom RBD simulation produced a non-finite final state.")
        end
    end

    return q_end, v_end
end

end # module
