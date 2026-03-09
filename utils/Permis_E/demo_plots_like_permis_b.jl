import Dionysos

const DI = Dionysos
const PR = DI.Problem

module NewSolver
include(joinpath(@__DIR__, "run_new_solver_on_permis_b.jl"))
end

include(joinpath(@__DIR__, "helpers", "helpers_plots.jl"))

function build_problem_and_params()
    system_cfg = NewSolver.PBRef.build_concrete_system()
    control_cfg = NewSolver.PBRef.build_control_problem()
    problem = PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
    return problem, system_cfg.params
end

function main(; outdir = joinpath(@__DIR__, "outputs"))
    mkpath(outdir)

    problem, params = build_problem_and_params()
    sim_cfg = NewSolver.PBRef.SimulationConfig()
    res = NewSolver.main()

    save_state_space_plots!(
        outdir,
        problem,
        res.candidate;
        cert_result = res.certification,
        show_ellipsoids = true,
        unwrap_angles = false,
        wrap_angles = true,
        periodic_dims = sim_cfg.periodic_dims,
        periodic_periods = sim_cfg.periodic_periods,
        periodic_start = sim_cfg.periodic_start,
    )

    gifpath = joinpath(outdir, "rollout.gif")
    plot_articulated_vehicle!(
        NewSolver.PBRef.AV,
        problem.system,
        params,
        res.candidate.x_traj,
        res.candidate.u_traj;
        giffile = gifpath,
        dt = res.candidate.Ts,
    )

    println("Saved: ", outdir)
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
