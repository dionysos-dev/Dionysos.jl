import Dionysos

const DI = Dionysos
const PR = DI.Problem
const ST = DI.System
const OP = DI.Optim

module PBRef
include(joinpath(@__DIR__, "..", "Permis_B", "permis_B.jl"))
end

include(joinpath(@__DIR__, "helpers", "helpers_plots.jl"))

function _build_problem(concrete_system)
    ctrl = PBRef.build_control_problem()
    return PR.OptimalControlProblem(
        concrete_system,
        ctrl.initial_set,
        ctrl.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

function _candidate_from_run(run_result)
    return OP.CandidateTrajectory(
        ST.Trajectory(run_result.state_list),
        ST.Trajectory(run_result.input_list);
        Ts = run_result.Δt,
        source = :permis_b_nominal,
        metadata = (;),
    )
end

function main(; outdir = joinpath(@__DIR__, "pdf"))
    mkpath(outdir)
    terminal_pdf = joinpath(outdir, "lmi_ellipsoide_terminale_12.pdf")
    isfile(terminal_pdf) && rm(terminal_pdf)

    sim_cfg = PBRef.SimulationConfig()
    lmi_cfg = PBRef.LMIConfig()

    run_result = PBRef.run_nominal_simulation(
        sim_cfg;
        save_plots = false,
        show_animation = false,
        use_cache = true,
        force_recompute = false,
    )
    run_result_lmi = PBRef.prepare_run_result_for_lmi(
        run_result;
        enable_unwrap = true,
        periodic_dims = sim_cfg.periodic_dims,
        periodic_periods = sim_cfg.periodic_periods,
    )

    problem = _build_problem(run_result.concrete_system)
    cand = _candidate_from_run(run_result)

    # plot12 + plot34 (same names as Permis_B output expectation)
    save_state_space_plots!(
        outdir,
        problem,
        cand;
        cert_result = nothing,
        show_ellipsoids = false,
        unwrap_angles = false,
        wrap_angles = true,
        periodic_dims = sim_cfg.periodic_dims,
        periodic_periods = sim_cfg.periodic_periods,
        periodic_start = sim_cfg.periodic_start,
        title12 = "state_space_12",
        title34 = "state_space_34",
    )

    # LMI diagnostics without terminal-ellipsoid PDF.
    diag = PBRef.run_lmi_diagnostic(
        run_result_lmi;
        cfg = lmi_cfg,
        output_dir = outdir,
        save_plots = false,
        unwrap_periodic_angles = false,
    )
    PBRef.plot_transitions_backward_chaine(
        run_result_lmi.state_list,
        diag.transitions.ellipsoides;
        output_dir = outdir,
        dims = (1, 2),
        filename = "lmi_transition_backward_12.pdf",
    )
    PBRef.plot_transitions_backward_chaine(
        run_result_lmi.state_list,
        diag.transitions.ellipsoides;
        output_dir = outdir,
        dims = (3, 4),
        title = "Test LMI - Chaine transitions backward (angles theta, phi)",
        filename = "lmi_transition_backward_34_angles.pdf",
    )

    # Optional overlay of certified ellipsoids on state-space plots.
    save_state_space_plots!(
        outdir,
        problem,
        cand;
        cert_result = diag,
        show_ellipsoids = true,
        unwrap_angles = false,
        wrap_angles = true,
        periodic_dims = sim_cfg.periodic_dims,
        periodic_periods = sim_cfg.periodic_periods,
        periodic_start = sim_cfg.periodic_start,
        title12 = "state_space_12 + ellipsoids",
        title34 = "state_space_34 + ellipsoids",
    )

    gifpath = joinpath(outdir, "centered_simu.gif")
    plot_articulated_vehicle!(
        PBRef.AV,
        run_result.concrete_system,
        run_result.params,
        cand.x_traj,
        cand.u_traj;
        giffile = gifpath,
        dt = 0.09,
        title = "centered_simu",
    )

    println("Saved files in: ", outdir)
    return (; outdir, run_result, diag)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
