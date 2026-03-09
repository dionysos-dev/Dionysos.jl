import Dionysos

const DI = Dionysos
const PR = DI.Problem
const OP = DI.Optim
const ST = DI.System
const UT = DI.Utils

module NewSolver
include(joinpath(@__DIR__, "run_new_solver_on_permis_b.jl"))
end

function build_problem()
    system_cfg = NewSolver.PBRef.build_concrete_system()
    control_cfg = NewSolver.PBRef.build_control_problem()
    return PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

function extract_tube(cert_result)
    cert_result === nothing && return nothing

    ells = UT.Ellipsoid[]

    if hasproperty(cert_result, :lmi_data) &&
       cert_result.lmi_data !== nothing &&
       hasproperty(cert_result.lmi_data, :ellipsoids)
        for E in cert_result.lmi_data.ellipsoids
            E isa UT.Ellipsoid && push!(ells, E)
        end
    end

    if isempty(ells) && hasproperty(cert_result, :steps)
        for rec in cert_result.steps
            hasproperty(rec, :ellipsoid) || continue
            hasproperty(rec, :status) && rec.status != :ok && continue
            E = rec.ellipsoid
            E === nothing && continue
            E isa UT.Ellipsoid && push!(ells, E)
        end
    end

    isempty(ells) && return nothing

    # Backward chain -> forward indexing approximation for x_k checks.
    return OP.TubeSpec(reverse(ells))
end

function main()
    problem = build_problem()
    res = NewSolver.main()

    us = collect(ST.enum_elems(res.candidate.u_traj))
    controller = OP.OpenLoopController(us)
    tube = extract_tube(res.certification)

    mc_cfg = OP.MonteCarloConfig(N = 200, seed = 1, keep_traces = 5)
    ro_cfg = OP.RolloutConfig(
        Ts = res.candidate.Ts,
        nstep = length(us),
        num_substeps = 5,
        keep_trace = false,
        stop_on_target = true,
        stop_on_violation = true,
    )

    mc_res = OP.monte_carlo(
        problem,
        controller;
        mc_cfg = mc_cfg,
        rollout_cfg = ro_cfg,
        tube = tube,
    )

    println("success_rate = ", mc_res.success_rate)
    println("violation_rate = ", mc_res.violation_rate)
    println("domain_violation_rate = ", mc_res.domain_violation_rate)
    println("tube_violation_rate = ", mc_res.tube_violation_rate)
    println("mean_time_to_target_s = ", mc_res.mean_time_to_target)

    return mc_res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
