# c'est fait à la va vite, il faut fixe ça demain dans la nouvelle version de dio
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import MathematicalSystems as MS
import Random

# Load the articulated vehicle benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the marche arriere benchmark.
"""
Base.@kwdef struct MarcheArriereConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 40

    terminal_radius::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.15, 0.15),
        IA.interval(-0.2, 0.2),
    )
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    verbose::Bool = false

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 5
    mppi_nsamples::Int = 1900
    mppi_niter::Int = 5 * 15
    mppi_λ::Float64 = 3.0
    mppi_noise_v::Float64 = 0.35
    mppi_noise_σ::Float64 = 0.15
end


######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(0.0, 0.0, -pi, -pi),
        SVector(12.0, 10.0, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(35.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(3.0, 4.0), SVector(9.0, 6.0)),
        UT.HyperRectangle(SVector(4.0, 3.0), SVector(8.0, 7.0)),
        UT.HyperRectangle(SVector(5.0, 2.0), SVector(7.0, 8.0)),
    ]

    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping()
    inputs_delta = [
        [-2.0, 0.0],
        [0.0, 0.0],
        [-1.0, 0.0],
        [-0.5, 0.0],
        [-1.5, 0.1], [-1.5, -0.1],
        [-1.0, 0.1], [-1.0, -0.1],
        [-0.5, 0.1], [-0.5, -0.1],
        [-1.5, 0.35], [-1.5, -0.35],
        [-1.0, 0.35], [-1.0, -0.35],
        [-0.5, 0.35], [-0.5, -0.35],
        [-1.5, 0.65], [-1.5, -0.65],
        [-1.0, 0.65], [-1.0, -0.65],
        [-0.5, 0.65], [-0.5, -0.65],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    x0 = SVector(6.0, 9.0, pi, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(5.25, 8.50, pi - deg2rad(4.0), -deg2rad(2.0)),
        SVector(6.45, 9.50, pi + deg2rad(4.0), deg2rad(2.0)),
    )

    target_set = UT.HyperRectangle(
        SVector(5.00, 0.50, -deg2rad(6.0), -deg2rad(4.0)),
        SVector(6.60, 1.5, deg2rad(6.0), deg2rad(4.0)),
    )

    return (; x0, initial_set, target_set)
end

function build_wrap_state(cfg::MarcheArriereConfig)
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

function periodic_state_error(x, xref, cfg::MarcheArriereConfig)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    return SVector{4, Float64}(e)
end

function project_input_to_domain(u, u_domain)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function build_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::MarcheArriereConfig;
    vehicle_module,
    input_mapping,
)
    seed_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        input_mapping,
        vehicle_module.jacobian_bound(system_cfg.params),
        periodicity_kwargs(cfg),
        cfg.nstep,
        _ -> control_cfg.x0;
        num_substeps = cfg.seed_num_substeps,
        trajectory_mode = cfg.seed_trajectory_mode,
    )

    seed_gen = OP.CenteredAbstractionGenerator(problem, seed_cfg)

    disc_system = ST.discretize_continuous_system(
        system_cfg.concrete_system,
        cfg.Δt;
        num_substeps = cfg.seed_num_substeps,
    )
    f_disc = MS.mapping(disc_system)
    wrap_state = build_wrap_state(cfg)
    target_center = UT.get_center(control_cfg.target_set)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    noise_sampler = function (rng, u, k)
        return [
            cfg.mppi_noise_v * Random.randn(rng),
            cfg.mppi_noise_σ * Random.randn(rng),
        ]
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)

    trajectory_cost = function (prob, cand)
        xs = collect(ST.enum_elems(cand.x_traj))
        us = collect(ST.enum_elems(cand.u_traj))

        BAD_COST = 1.0e12
        J = 0.0
        w_step = 2.0
        w_pos = 1.0
        w_ang = 0.15
        w_u = 0.05
        w_du = 0.5
        w_ddu = 0.1

        hit_target = false
        hit_index = length(xs)

        for k in eachindex(xs)
            xw = wrap_state(xs[k])

            if !(xw ∈ prob.system.X)
                return BAD_COST
            end

            if xw ∈ prob.target_set
                hit_target = true
                hit_index = k
                break
            end

            e = periodic_state_error(xw, target_center, cfg)
            J += w_step
            J += w_pos * (e[1]^2 + e[2]^2)
            J += w_ang * (e[3]^2 + e[4]^2)
        end

        if !hit_target
            xT = wrap_state(last(xs))
            eT = periodic_state_error(xT, target_center, cfg)
            J += 500.0 * (eT[1]^2 + eT[2]^2)
            J += 80.0 * (eT[3]^2 + eT[4]^2)
            J += 1.0e4
        end

        last_u_index = min(length(us), max(hit_index - 1, 0))

        for k in 1:last_u_index
            u = us[k]
            J += w_u * sum(abs2, u)

            if k >= 2
                du = u .- us[k - 1]
                J += w_du * sum(abs2, du)
            end

            if k >= 3
                ddu = u .- 2 .* us[k - 1] .+ us[k - 2]
                J += w_ddu * sum(abs2, ddu)
            end
        end

        return J
    end

    success_fun = function (prob, cand)
        xT = wrap_state(last(ST.enum_elems(cand.x_traj)))
        return xT ∈ prob.target_set
    end

    mppi_cfg = OP.MPPIConfig(
        cfg.Δt,
        cfg.nstep,
        cfg.mppi_nsamples,
        cfg.mppi_niter,
        cfg.mppi_λ,
        _ -> control_cfg.x0,
        discrete_dynamics,
        noise_sampler,
        project_input,
        trajectory_cost,
        success_fun,
        (prob, x) -> wrap_state(x),
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function save_named_state_space_plots!(
    plots_dir::AbstractString,
    basename::AbstractString,
    problem,
    candidate;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims,
    periodic_periods,
    periodic_start,
    title12::AbstractString,
    title34::AbstractString,
)
    mktempdir() do tmp_dir
        save_state_space_plots!(
            tmp_dir,
            problem,
            candidate;
            cert_result = cert_result,
            show_ellipsoids = show_ellipsoids,
            unwrap_angles = unwrap_angles,
            wrap_angles = wrap_angles,
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
            periodic_start = periodic_start,
            title12 = title12,
            title34 = title34,
        )

        cp(
            joinpath(tmp_dir, "state_space_12.pdf"),
            joinpath(plots_dir, "$(basename)_12.pdf");
            force = true,
        )
        cp(
            joinpath(tmp_dir, "state_space_34.pdf"),
            joinpath(plots_dir, "$(basename)_34.pdf");
            force = true,
        )
    end

    return nothing
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    run_result = run_vehicle_benchmark(
        cfg;
        scenario_name = "cercle_mppi",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
        generator_builder = build_mppi_generator,
    )

    warmup_candidate = OP.get_seed(run_result.solver.generator)
    warmup_candidate === nothing && error("MPPI warmup trajectory is missing.")

    plots_dir = run_result.outputs.plots_dir
    problem = run_result.problem
    cert_result = run_result.result.certification

    save_named_state_space_plots!(
        plots_dir,
        "abstract_traj_used_as_warmup",
        problem,
        warmup_candidate;
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "Abstract traj used as warmup",
        title34 = "Abstract traj used as warmup",
    )

    save_named_state_space_plots!(
        plots_dir,
        "mppi_candidate_traj",
        problem,
        run_result.nominal_candidate;
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "MPPI candidate traj",
        title34 = "MPPI candidate traj",
    )

    save_named_state_space_plots!(
        plots_dir,
        "ellipsoids",
        problem,
        run_result.nominal_candidate;
        cert_result = cert_result,
        show_ellipsoids = true,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "cercle_mppi (x,y)",
        title34 = "cercle_mppi (theta,phi)",
    )

    return run_result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
