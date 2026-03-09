import Dionysos
import StaticArrays: SVector

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PR = DI.Problem
const OP = DI.Optim

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

Base.@kwdef struct NominalConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 700
end

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi),
        SVector(10.0, 9.0, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(50.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 6.0), SVector(10.0, 9.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)
    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_domain()
    inputs_delta = [
        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return DO.CustomList(inputs)
end

function build_control_problem()
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    initial_set = UT.HyperRectangle(
        SVector(-1.0, -1.0, -0.4, -0.4),
        SVector(1.0, 1.0, 0.4, 0.4),
    )

    target_set = UT.HyperRectangle(
        SVector(9.0, 5.0, pi - 5 * (pi / 180), -5 * (pi / 180)),
        SVector(10.0, 6.0, pi + 5 * (pi / 180), 5 * (pi / 180)),
    )

    return (; x0, initial_set, target_set)
end

periodicity_kwargs(cfg::NominalConfig) = (;
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

function get_nominal_traje(cfg::NominalConfig = NominalConfig())
    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    Udom = build_input_domain()

    problem = PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )

    gen_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        Udom,
        AV.jacobian_bound(system_cfg.params),
        periodicity_kwargs(cfg),
        cfg.nstep,
        _ -> control_cfg.x0,
    )

    gen = OP.CenteredAbstractionGenerator{typeof(problem), typeof(gen_cfg), Any, Any}(
        problem,
        gen_cfg,
        nothing,
        nothing,
        false,
        0.0,
    )

    OP.set_problem!(gen, problem)
    OP.generate!(gen)

    cand = OP.get_trajectory(gen)
    return (;
        success = OP.get_success(gen),
        solve_time_sec = OP.get_solve_time(gen),
        candidate = cand,
        x_traj = cand === nothing ? nothing : cand.x_traj,
        u_traj = cand === nothing ? nothing : cand.u_traj,
        problem = problem,
        generator = gen,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    res = get_nominal_traje()
    println("success = ", res.success)
    println("solve_time_sec = ", res.solve_time_sec)
    println("candidate = ", res.candidate)
    println("is CandidateTrajectory = ", res.candidate isa OP.CandidateTrajectory)
    if res.candidate !== nothing
        println("horizon = ", OP.horizon(res.candidate))
        println("n_states = ", OP.n_states(res.candidate))
    end
end
