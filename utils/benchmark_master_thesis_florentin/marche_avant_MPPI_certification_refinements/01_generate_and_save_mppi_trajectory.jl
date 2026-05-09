include(joinpath(@__DIR__, "run_marche_avant.jl"))

using JLD2
using StaticArrays: SVector

const RESULTS_DIR = joinpath(@__DIR__, "results")
const TRAJECTORY_PATH = joinpath(RESULTS_DIR, "saved_mppi_trajectory.jld2")

function _trajectory_matrix(traj)
    elems = collect(ST.enum_elems(traj))
    return reduce(hcat, collect.(Float64, elems))
end

function main()
    mkpath(RESULTS_DIR)

    cfg = MarcheAvantConfig(;
        output_root = joinpath(RESULTS_DIR, "generation_outputs"),
        plot_gif = false,
    )

    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)

    gen = build_mppi_generator(
        problem,
        system_cfg,
        control_cfg,
        cfg;
        vehicle_module = AV,
        input_mapping = build_input_mapping(),
    )

    OP.set_problem!(gen, problem)
    OP.generate!(gen)

    candidate = OP.get_trajectory(gen)
    candidate === nothing && error("MPPI did not produce a candidate trajectory.")

    x_traj = _trajectory_matrix(candidate.x_traj)
    u_traj = _trajectory_matrix(candidate.u_traj)
    Ts = candidate.Ts
    source = candidate.source
    generator_success = OP.get_success(gen)
    generator_solve_time = OP.get_solve_time(gen)
    generator_diagnostics = gen.diagnostics
    config_metadata = (;
        Δt = cfg.Δt,
        nstep = cfg.nstep,
        terminal_radius = cfg.terminal_radius,
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        mppi_nsamples = cfg.mppi_nsamples,
        mppi_niter = cfg.mppi_niter,
        mppi_λ = cfg.mppi_λ,
        mppi_noise_v = cfg.mppi_noise_v,
        mppi_noise_σ = cfg.mppi_noise_σ,
    )

    jldsave(
        TRAJECTORY_PATH;
        x_traj,
        u_traj,
        Ts,
        source,
        generator_success,
        generator_solve_time,
        generator_diagnostics,
        config_metadata,
    )

    println("saved MPPI trajectory: ", TRAJECTORY_PATH)
    println("generator_success = ", generator_success)
    println("horizon = ", size(u_traj, 2))
    return TRAJECTORY_PATH
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
