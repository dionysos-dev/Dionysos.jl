include(joinpath(@__DIR__, "04_final_scaling_validation_sweep.jl"))

using DataFrames
using LinearAlgebra
using Plots
using Random
using Statistics
using StaticArrays: SVector

const TRANSITION_AUDIT_DIR = joinpath(RESULTS_DIR, "transition_audit")
const TRANSITION_AUDIT_PLOTS_DIR = joinpath(TRANSITION_AUDIT_DIR, "plots")
const TRANSITION_AUDIT_SUMMARY = joinpath(TRANSITION_AUDIT_DIR, "transition_audit_summary.csv")
const TRANSITION_AUDIT_PER_TRANSITION =
    joinpath(TRANSITION_AUDIT_DIR, "transition_audit_per_transition.csv")
const TRANSITION_AUDIT_WORST_SAMPLES =
    joinpath(TRANSITION_AUDIT_DIR, "transition_audit_worst_samples.csv")
const TRANSITION_AUDIT_GLOBAL_REASONS =
    joinpath(TRANSITION_AUDIT_DIR, "global_stress_failure_reasons.csv")
const TRANSITION_AUDIT_SIMULATORS =
    joinpath(TRANSITION_AUDIT_DIR, "simulator_comparison.csv")
const TRANSITION_AUDIT_DIAGNOSTICS =
    joinpath(TRANSITION_AUDIT_DIR, "transition_audit_diagnostics.md")

const audit_rng_seed = 20250428
const audit_uniform_samples =
    parse(Int, get(ENV, "TRANSITION_AUDIT_UNIFORM_SAMPLES", "10000"))
const audit_boundary_samples =
    parse(Int, get(ENV, "TRANSITION_AUDIT_BOUNDARY_SAMPLES", "2000"))
const audit_global_samples =
    parse(Int, get(ENV, "TRANSITION_AUDIT_GLOBAL_SAMPLES", "1000"))
const audit_tol = 1.0e-7
const audit_tols = (1.0e-7, 1.0e-6, 1.0e-5)

_audit_vec_string(x) = x === nothing ? "nothing" : join(string.(Float64.(x)), ";")

function _audit_interval_radii(Δ)
    return [max(abs(Float64(IA.inf(I))), abs(Float64(IA.sup(I)))) for I in collect(Δ)]
end

function _audit_bounds_string(set)
    if hasproperty(set, :lb) && hasproperty(set, :ub)
        return join(string.(Float64.(set.lb)), ";") * "|" * join(string.(Float64.(set.ub)), ";")
    end
    return string(typeof(set))
end

function _audit_kappa_matrix(kappa)
    hasproperty(kappa, :A) || error("Unsupported kappa type without A field: $(typeof(kappa))")
    K = Matrix{Float64}(getproperty(kappa, :A))
    b =
        hasproperty(kappa, :b) ? vec(Float64.(getproperty(kappa, :b))) :
        (hasproperty(kappa, :c) ? vec(Float64.(getproperty(kappa, :c))) :
         error("Unsupported kappa type without b/c offset: $(typeof(kappa))"))
    return K, b
end

function _audit_ellipsoid_value(E, x)
    dx = collect(Float64, x) .- collect(Float64, E.c)
    return Float64(dx' * Matrix(E.P) * dx)
end

function _audit_periodic_error(x, y, cfg)
    e = collect(Float64, x) .- collect(Float64, y)
    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end
    return e
end

function _audit_sample_ellipsoid(E; n_uniform, n_boundary, rng)
    P = Symmetric(Matrix{Float64}(E.P))
    U = cholesky(P).U
    c = collect(Float64, E.c)
    nx = length(c)
    points = Vector{Vector{Float64}}()
    labels = String[]

    push!(points, copy(c))
    push!(labels, "center")

    for i in 1:n_uniform
        v = randn(rng, nx)
        while norm(v) <= 1.0e-14
            v = randn(rng, nx)
        end
        z = (rand(rng)^(1 / nx) / norm(v)) .* v
        push!(points, collect(c .+ U \ z))
        push!(labels, "uniform")
    end

    for i in 1:n_boundary
        v = randn(rng, nx)
        while norm(v) <= 1.0e-14
            v = randn(rng, nx)
        end
        z = v ./ norm(v)
        push!(points, collect(c .+ U \ z))
        push!(labels, "boundary")
    end

    Q = Symmetric(inv(P))
    for i in 1:nx
        e = zeros(nx)
        e[i] = 1.0
        d = Q * e
        denom = sqrt(e' * Q * e)
        denom <= 1.0e-14 && continue
        push!(points, collect(c .+ d ./ denom))
        push!(labels, "coord_extreme_plus_$i")
        push!(points, collect(c .- d ./ denom))
        push!(labels, "coord_extreme_minus_$i")
    end

    return points, labels
end

function _audit_sample_uniform_ellipsoid(E, n; rng)
    P = Symmetric(Matrix{Float64}(E.P))
    U = cholesky(P).U
    c = collect(Float64, E.c)
    nx = length(c)
    points = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        v = randn(rng, nx)
        while norm(v) <= 1.0e-14
            v = randn(rng, nx)
        end
        z = (rand(rng)^(1 / nx) / norm(v)) .* v
        points[i] = collect(c .+ U \ z)
    end
    return points
end

function _audit_box_inclusion(E, xbar, box_radii; tol = audit_tol)
    P = Symmetric(Matrix{Float64}(E.P))
    Q = Symmetric(inv(P))
    ell_radii = sqrt.(max.(diag(Q), 0.0))
    center_distance = abs.(collect(Float64, E.c) .- collect(Float64, xbar))
    margins = center_distance .+ ell_radii .- box_radii
    worst = argmax(margins)
    return (;
        pass = maximum(margins) <= tol,
        max_margin = maximum(margins),
        worst_dim = worst,
        center_distance,
        ellipsoid_radius = ell_radii,
        box_radius = box_radii,
    )
end

function _audit_input_inclusion(E, kappa, U; tol = audit_tol)
    K, b = _audit_kappa_matrix(kappa)
    P = Symmetric(Matrix{Float64}(E.P))
    Q = Symmetric(inv(P))
    u_center = K * collect(Float64, E.c) + b
    input_radius = [sqrt(max(Float64(K[j, :]' * Q * K[j, :]), 0.0)) for j in 1:size(K, 1)]

    if hasproperty(U, :lb) && hasproperty(U, :ub)
        lb = collect(Float64, U.lb)
        ub = collect(Float64, U.ub)
        U_center = (lb .+ ub) ./ 2
        U_radius = (ub .- lb) ./ 2
        margins = abs.(u_center .- U_center) .+ input_radius .- U_radius
        worst = argmax(margins)
        return (;
            pass = maximum(margins) <= tol,
            max_margin = maximum(margins),
            worst_dim = worst,
            nominal_u = u_center,
            input_radius,
            bounds = _audit_bounds_string(U),
        )
    end

    return (;
        pass = false,
        max_margin = Inf,
        worst_dim = missing,
        nominal_u = u_center,
        input_radius,
        bounds = string(typeof(U)),
    )
end

function _audit_independent_rk4_map(system_cfg, Ts, num_substeps)
    f = AV.dynamic(system_cfg.params)
    h = Float64(Ts) / num_substeps
    return function (x, u)
        xv = SVector{4, Float64}(x)
        uv = SVector{2, Float64}(u)
        for _ in 1:num_substeps
            k1 = f(xv, uv)
            k2 = f(xv + 0.5 * h * k1, uv)
            k3 = f(xv + 0.5 * h * k2, uv)
            k4 = f(xv + h * k3, uv)
            xv = xv + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        end
        return collect(xv)
    end
end

function _audit_nominal_indexing_rows(result, candidate, cfg)
    chain = _build_certified_kappa_chain(result)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    rows = NamedTuple[]
    for k in chain.k_sequence
        E = chain.ellipsoid_by_k[k]
        Enext = chain.ellipsoid_by_k[k + 1]
        K, b = _audit_kappa_matrix(chain.kappa_by_k[k])
        push!(
            rows,
            (;
                k,
                source_center_minus_nominal_norm = norm(_audit_periodic_error(E.c, xs[k], cfg)),
                target_center_minus_nominal_norm = norm(_audit_periodic_error(Enext.c, xs[k + 1], cfg)),
                kappa_at_source_center_minus_u_norm = norm(K * collect(E.c) + b - collect(us[k])),
                source_center = _audit_vec_string(E.c),
                source_nominal = _audit_vec_string(xs[k]),
                target_center = _audit_vec_string(Enext.c),
                target_nominal = _audit_vec_string(xs[k + 1]),
                nominal_control = _audit_vec_string(us[k]),
            ),
        )
    end
    return rows
end

function _audit_transition_empirical(
    config_name,
    run,
    problem,
    system_cfg,
    candidate,
    cfg;
    rng,
)
    chain = _build_certified_kappa_chain(run.result)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    maps = _build_periodic_maps(
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
    )
    fA = _build_discrete_system_map(problem.system, candidate.Ts; num_substeps = cfg.symbolic_rk4_substeps)
    fB = _audit_independent_rk4_map(system_cfg, candidate.Ts, cfg.symbolic_rk4_substeps)
    box_radii = _audit_interval_radii(cfg.ΔX)
    U = problem.system.U

    per_transition = NamedTuple[]
    worst_rows = NamedTuple[]
    simulator_rows = NamedTuple[]

    for k in chain.k_sequence
        E = chain.ellipsoid_by_k[k]
        Enext = chain.ellipsoid_by_k[k + 1]
        kappa = chain.kappa_by_k[k]
        points, labels = _audit_sample_ellipsoid(
            E;
            n_uniform = audit_uniform_samples,
            n_boundary = audit_boundary_samples,
            rng,
        )
        box = _audit_box_inclusion(E, xs[k], box_radii)
        input = _audit_input_inclusion(E, kappa, U)

        source_fail = 0
        box_fail = 0
        input_fail = 0
        next_fail = Dict(t => 0 for t in audit_tols)
        values_next = Float64[]
        sim_diffs = Float64[]
        worst_value = -Inf
        worst_x = Float64[]
        worst_u = Float64[]
        worst_xnext = Float64[]
        worst_label = ""
        max_input_violation_sample = 0.0

        for (xraw, label) in zip(points, labels)
            x = maps.lift_state_near_reference(maps.wrap_state(xraw), xs[k])
            vsrc = _audit_ellipsoid_value(E, x)
            vsrc > 1 + audit_tol && (source_fail += 1)

            center_dist = abs.(collect(x) .- collect(xs[k]))
            any(center_dist .> box_radii .+ audit_tol) && (box_fail += 1)

            u = _eval_kappa_controller(kappa, x)
            if hasproperty(U, :lb) && hasproperty(U, :ub)
                viol = maximum(max.(collect(U.lb) .- u, u .- collect(U.ub), 0.0))
                max_input_violation_sample = max(max_input_violation_sample, viol)
                viol > audit_tol && (input_fail += 1)
            end

            xnextA_phys = maps.wrap_state(fA(maps.wrap_state(x), u))
            xnextB_phys = maps.wrap_state(fB(maps.wrap_state(x), u))
            diffAB = norm(_audit_periodic_error(xnextA_phys, xnextB_phys, cfg))
            push!(sim_diffs, diffAB)

            xnext = maps.lift_state_near_reference(xnextA_phys, xs[k + 1])
            vnext = _audit_ellipsoid_value(Enext, xnext)
            push!(values_next, vnext)
            for tol in audit_tols
                vnext > 1 + tol && (next_fail[tol] += 1)
            end
            if vnext > worst_value
                worst_value = vnext
                worst_x = collect(x)
                worst_u = collect(u)
                worst_xnext = collect(xnext)
                worst_label = label
            end
        end

        n = length(points)
        p99 = quantile(values_next, 0.99)
        push!(
            per_transition,
            (;
                configuration = config_name,
                k,
                lmi_step_status = string(run.result.steps[k].status),
                box_inclusion_pass = box.pass,
                max_box_margin = box.max_margin,
                worst_box_dimension = box.worst_dim,
                center_distance_per_dim = _audit_vec_string(box.center_distance),
                ellipsoid_radius_per_dim = _audit_vec_string(box.ellipsoid_radius),
                box_radius_per_dim = _audit_vec_string(box.box_radius),
                input_inclusion_pass = input.pass,
                max_input_margin = input.max_margin,
                worst_input_dimension = input.worst_dim,
                nominal_u = _audit_vec_string(input.nominal_u),
                input_radius_from_feedback = _audit_vec_string(input.input_radius),
                U_bounds = input.bounds,
                n_samples = n,
                source_membership_failures = source_fail,
                box_membership_failures = box_fail,
                input_membership_failures = input_fail,
                next_ellipsoid_failures_tol_1e_7 = next_fail[1.0e-7],
                next_ellipsoid_failures_tol_1e_6 = next_fail[1.0e-6],
                next_ellipsoid_failures_tol_1e_5 = next_fail[1.0e-5],
                max_next_ellipsoid_value = maximum(values_next),
                p99_next_ellipsoid_value = p99,
                mean_next_ellipsoid_value = mean(values_next),
                simulator_max_difference = maximum(sim_diffs),
                simulator_mean_difference = mean(sim_diffs),
                simulator_p99_difference = quantile(sim_diffs, 0.99),
                nominal_control = _audit_vec_string(us[k]),
            ),
        )
        push!(
            worst_rows,
            (;
                configuration = config_name,
                k,
                worst_sample_label = worst_label,
                worst_next_ellipsoid_value = worst_value,
                worst_sample_x = _audit_vec_string(worst_x),
                worst_sample_u = _audit_vec_string(worst_u),
                worst_sample_x_next = _audit_vec_string(worst_xnext),
            ),
        )
        push!(
            simulator_rows,
            (;
                configuration = config_name,
                k,
                n_samples = n,
                max_difference = maximum(sim_diffs),
                mean_difference = mean(sim_diffs),
                p99_difference = quantile(sim_diffs, 0.99),
                time_step = candidate.Ts,
                num_substeps = cfg.symbolic_rk4_substeps,
                note = "Simulator A uses ST.discretize_continuous_system; simulator B is local RK4 over AV.dynamic.",
            ),
        )
    end
    return (; per_transition, worst_rows, simulator_rows)
end

function _audit_global_stress(
    config_name,
    run,
    problem,
    system_cfg,
    candidate,
    cfg;
    rng,
)
    chain = _build_certified_kappa_chain(run.result)
    xs = collect(ST.enum_elems(candidate.x_traj))
    maps = _build_periodic_maps(
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
    )
    fA = _build_discrete_system_map(problem.system, candidate.Ts; num_substeps = cfg.symbolic_rk4_substeps)
    fB = _audit_independent_rk4_map(system_cfg, candidate.Ts, cfg.symbolic_rk4_substeps)
    box_radii = _audit_interval_radii(cfg.ΔX)
    U = problem.system.U
    domain_set = _domain_set(problem.system)
    obstacle_set = _obstacle_set(problem.system)
    samples = _audit_sample_uniform_ellipsoid(chain.initial_ellipsoid, audit_global_samples; rng)
    rows = NamedTuple[]

    for (sample_id, x0) in enumerate(samples)
        status = "success"
        reason = "success"
        first_failure_step = missing
        max_next_value = -Inf
        max_sim_diff = 0.0
        x = maps.wrap_state(x0)

        x_chain0 = maps.lift_state_near_reference(x, xs[first(chain.k_sequence)])
        if _audit_ellipsoid_value(chain.initial_ellipsoid, x_chain0) > 1 + audit_tol
            status = "failure"
            reason = "initial_not_in_E0"
            first_failure_step = first(chain.k_sequence)
        end

        for k in chain.k_sequence
            status == "success" || break
            E = chain.ellipsoid_by_k[k]
            Enext = chain.ellipsoid_by_k[k + 1]
            x_chain = maps.lift_state_near_reference(maps.wrap_state(x), xs[k])

            if _audit_ellipsoid_value(E, x_chain) > 1 + audit_tol
                status = "failure"
                reason = "left_source_ellipsoid_before_step"
                first_failure_step = k
                break
            end
            if any(abs.(collect(x_chain) .- collect(xs[k])) .> box_radii .+ audit_tol)
                status = "failure"
                reason = "left_Xbar"
                first_failure_step = k
                break
            end

            u = _eval_kappa_controller(chain.kappa_by_k[k], x_chain)
            if hasproperty(U, :lb) && hasproperty(U, :ub)
                viol = maximum(max.(collect(U.lb) .- u, u .- collect(U.ub), 0.0))
                if viol > audit_tol
                    status = "failure"
                    reason = "input_outside_U"
                    first_failure_step = k
                    break
                end
            end

            xnextA = maps.wrap_state(fA(maps.wrap_state(x_chain), u))
            xnextB = maps.wrap_state(fB(maps.wrap_state(x_chain), u))
            sim_diff = norm(_audit_periodic_error(xnextA, xnextB, cfg))
            max_sim_diff = max(max_sim_diff, sim_diff)
            if sim_diff > 1.0e-8
                status = "failure"
                reason = "simulator_mismatch"
                first_failure_step = k
                break
            end

            xnext_chain = maps.lift_state_near_reference(xnextA, xs[k + 1])
            vnext = _audit_ellipsoid_value(Enext, xnext_chain)
            max_next_value = max(max_next_value, vnext)
            if vnext > 1 + audit_tol
                status = "failure"
                reason = "next_not_in_target_ellipsoid"
                first_failure_step = k + 1
                break
            end
            if _point_in_any_set(xnextA, obstacle_set)
                status = "failure"
                reason = "obstacle_violation"
                first_failure_step = k + 1
                break
            end
            if domain_set !== nothing && !(xnextA ∈ domain_set)
                status = "failure"
                reason = "domain_violation"
                first_failure_step = k + 1
                break
            end
            x = xnextA
        end

        if status == "success" && !(maps.wrap_state(x) ∈ problem.target_set)
            status = "failure"
            reason = "target_not_reached"
            first_failure_step = last(chain.k_sequence) + 1
        end

        push!(
            rows,
            (;
                configuration = config_name,
                sample_id,
                status,
                first_failure_reason = reason,
                first_failure_step,
                max_next_ellipsoid_value = max_next_value,
                max_simulator_difference = max_sim_diff,
            ),
        )
    end
    return rows
end

function _audit_summary_for_config(
    config_name,
    run,
    per_rows,
    global_rows,
    simulator_rows;
    box_multiplier_X,
    state_scaling_name,
    state_scaling,
)
    lmi_feasible = run.result.success
    all_box = all(r.box_inclusion_pass for r in per_rows)
    all_input = all(r.input_inclusion_pass for r in per_rows)
    empirical_pass = all(r.next_ellipsoid_failures_tol_1e_7 == 0 for r in per_rows)
    simulators_ok = all(r.max_difference <= 1.0e-8 for r in simulator_rows)
    global_success_rate = count(r -> r.status == "success", global_rows) / length(global_rows)
    global_ok = global_success_rate >= 0.99

    verdict =
        !all_box ? "NOT FORMAL" :
        !all_input ? "NOT FORMAL" :
        !empirical_pass ? "BUG SUSPECTED" :
        !global_ok ? "GLOBAL STRESS TEST BUG SUSPECTED" :
        !simulators_ok ? "BUG SUSPECTED" : "FORMAL"

    reason_counts = join(
        ["$(r):$(count(row -> row.first_failure_reason == r, global_rows))" for r in sort(unique([row.first_failure_reason for row in global_rows]))],
        ";",
    )
    hist = join(
        [
            "$(s):$(count(row -> !ismissing(row.first_failure_step) && row.first_failure_step == s, global_rows))" for
            s in sort(unique(collect(skipmissing([row.first_failure_step for row in global_rows]))))
        ],
        ";",
    )

    return (;
        configuration = config_name,
        state_scaling_name,
        state_scaling_vector = _audit_vec_string(state_scaling),
        box_multiplier_X,
        lmi_feasible,
        certified_transitions = count(step -> step.status == :ok, run.result.steps),
        horizon = OP.horizon(candidate_global_ref[]),
        box_inclusion_all_transitions = all_box,
        input_inclusion_all_transitions = all_input,
        transition_empirical_test_pass = empirical_pass,
        simulator_A_B_consistent = simulators_ok,
        global_stress_test_consistent = global_ok,
        global_success_rate,
        max_box_margin = maximum(r.max_box_margin for r in per_rows),
        max_input_margin = maximum(r.max_input_margin for r in per_rows),
        max_next_ellipsoid_value = maximum(r.max_next_ellipsoid_value for r in per_rows),
        total_next_failures_tol_1e_7 = sum(r.next_ellipsoid_failures_tol_1e_7 for r in per_rows),
        total_box_membership_failures = sum(r.box_membership_failures for r in per_rows),
        total_input_membership_failures = sum(r.input_membership_failures for r in per_rows),
        simulator_max_difference = maximum(r.max_difference for r in simulator_rows),
        failure_reason_counts = reason_counts,
        first_failure_step_histogram = hist,
        final_verdict = verdict,
    )
end

# Tiny mutable holder to avoid threading `candidate` through the summary row only.
const candidate_global_ref = Ref{Any}(nothing)

function _audit_write_plots(per_transition_rows, global_rows, simulator_rows)
    mkpath(TRANSITION_AUDIT_PLOTS_DIR)
    per = DataFrame(per_transition_rows)
    glob = DataFrame(global_rows)
    sim = DataFrame(simulator_rows)

    p = plot(; xlabel = "transition k", ylabel = "max next ellipsoid value", yscale = :log10, title = "Transition audit: max v_{k+1}")
    for cfg in unique(per.configuration)
        df = per[per.configuration .== cfg, :]
        plot!(p, df.k, df.max_next_ellipsoid_value; label = cfg, marker = :circle)
    end
    savefig(p, joinpath(TRANSITION_AUDIT_PLOTS_DIR, "max_next_ellipsoid_value_by_transition.png"))

    p = plot(; xlabel = "transition k", ylabel = "violations", title = "Transition audit: empirical violations")
    for cfg in unique(per.configuration)
        df = per[per.configuration .== cfg, :]
        plot!(p, df.k, df.next_ellipsoid_failures_tol_1e_7; label = cfg, marker = :circle)
    end
    savefig(p, joinpath(TRANSITION_AUDIT_PLOTS_DIR, "violations_by_transition.png"))

    failed = glob[glob.status .!= "success", :]
    if !isempty(failed)
        p = histogram(
            skipmissing(failed.first_failure_step);
            bins = 1:34,
            xlabel = "first failure step",
            ylabel = "count",
            title = "Global stress first failure step",
            label = "",
        )
        savefig(p, joinpath(TRANSITION_AUDIT_PLOTS_DIR, "global_first_failure_step_histogram.png"))
    end

    p = plot(; xlabel = "transition k", ylabel = "Simulator A/B max difference", yscale = :log10, title = "Simulator comparison")
    for cfg in unique(sim.configuration)
        df = sim[sim.configuration .== cfg, :]
        plot!(p, df.k, df.max_difference .+ eps(); label = cfg, marker = :circle)
    end
    savefig(p, joinpath(TRANSITION_AUDIT_PLOTS_DIR, "simulator_A_vs_B_difference.png"))
end

function _audit_write_diagnostics(summary_rows, indexing_rows)
    df = DataFrame(summary_rows)
    open(TRANSITION_AUDIT_DIAGNOSTICS, "w") do io
        println(io, "# Transition-level certificate audit\n")
        println(io, "Monte Carlo stress tests are empirical diagnostics, not formal proofs. A transition is treated as locally formal only when the LMI solved, the source ellipsoid is inside the local linearization box, the affine feedback maps the full source ellipsoid inside U, and the direct transition-level replay has no violations at the selected tolerance.\n")
        println(io, "## Indexing diagnostics")
        for row in indexing_rows[1:min(end, 5)]
            println(io, "- k=$(row.k): ||E_k.c - x_k||=$(row.source_center_minus_nominal_norm), ||E_{k+1}.c - x_{k+1}||=$(row.target_center_minus_nominal_norm), ||kappa(E_k.c)-u_k||=$(row.kappa_at_source_center_minus_u_norm)")
        end
        println(io, "\n## Verdicts")
        for row in eachrow(df)
            println(io, "- $(row.configuration): $(row.final_verdict), global_success=$(row.global_success_rate), box_all=$(row.box_inclusion_all_transitions), input_all=$(row.input_inclusion_all_transitions), empirical_pass=$(row.transition_empirical_test_pass), sim_ok=$(row.simulator_A_B_consistent)")
        end
    end
end

function main()
    mkpath(TRANSITION_AUDIT_DIR)
    mkpath(TRANSITION_AUDIT_PLOTS_DIR)
    isfile(TRAJECTORY_PATH) || error("Missing saved MPPI trajectory. Run 01_generate_and_save_mppi_trajectory.jl first.")

    cfg0 = MarcheAvantConfig(; output_root = RESULTS_DIR, plot_gif = false)
    candidate, data = _load_candidate()
    candidate_global_ref[] = candidate
    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)
    scaling_pairs = Dict(String(k) => v for (k, v) in _scaling_candidates(data["x_traj"], cfg0))

    selected = [
        ("none", 1.0, scaling_pairs["none"]),
        ("none", 2.0, scaling_pairs["none"]),
        ("none", 3.0, scaling_pairs["none"]),
        ("trajectory_std", 1.0, scaling_pairs["trajectory_std"]),
        ("trajectory_std", 2.0, scaling_pairs["trajectory_std"]),
        ("trajectory_std", 3.0, scaling_pairs["trajectory_std"]),
        ("trajectory_range", 1.0, scaling_pairs["trajectory_range"]),
        ("trajectory_range", 2.0, scaling_pairs["trajectory_range"]),
        ("trajectory_range", 3.0, scaling_pairs["trajectory_range"]),
    ]

    summary_rows = NamedTuple[]
    per_transition_rows = NamedTuple[]
    worst_rows = NamedTuple[]
    simulator_rows = NamedTuple[]
    global_rows = NamedTuple[]
    indexing_rows_all = NamedTuple[]

    println("Transition audit settings: uniform=$(audit_uniform_samples), boundary=$(audit_boundary_samples), global=$(audit_global_samples)")
    for (name, mx, scaling) in selected
        cfg = _cfg_with_boxes(cfg0; mx = mx, mu = 1.0)
        config_name = "$(name)__box_$(mx)"
        println("[audit] certifying $(config_name)")
        run = _certify_named(Symbol(name), scaling, problem, system_cfg, cfg, candidate)
        if !run.result.success
            push!(
                summary_rows,
                (;
                    configuration = config_name,
                    state_scaling_name = name,
                    state_scaling_vector = _audit_vec_string(scaling),
                    box_multiplier_X = mx,
                    lmi_feasible = false,
                    certified_transitions = count(step -> step.status == :ok, run.result.steps),
                    horizon = OP.horizon(candidate),
                    box_inclusion_all_transitions = false,
                    input_inclusion_all_transitions = false,
                    transition_empirical_test_pass = false,
                    simulator_A_B_consistent = missing,
                    global_stress_test_consistent = false,
                    global_success_rate = 0.0,
                    max_box_margin = missing,
                    max_input_margin = missing,
                    max_next_ellipsoid_value = missing,
                    total_next_failures_tol_1e_7 = missing,
                    total_box_membership_failures = missing,
                    total_input_membership_failures = missing,
                    simulator_max_difference = missing,
                    failure_reason_counts = "certification_failed",
                    first_failure_step_histogram = "",
                    final_verdict = "NOT FORMAL",
                ),
            )
            continue
        end

        append!(indexing_rows_all, [(; configuration = config_name, row...) for row in _audit_nominal_indexing_rows(run.result, candidate, cfg)])

        rng = Random.MersenneTwister(audit_rng_seed + round(Int, 100 * mx) + hash(name) % 10000)
        audit = _audit_transition_empirical(config_name, run, problem, system_cfg, candidate, cfg; rng)
        append!(per_transition_rows, audit.per_transition)
        append!(worst_rows, audit.worst_rows)
        append!(simulator_rows, audit.simulator_rows)

        global_audit = _audit_global_stress(config_name, run, problem, system_cfg, candidate, cfg; rng)
        append!(global_rows, global_audit)
        push!(
            summary_rows,
            _audit_summary_for_config(
                config_name,
                run,
                audit.per_transition,
                global_audit,
                audit.simulator_rows;
                box_multiplier_X = mx,
                state_scaling_name = name,
                state_scaling = scaling,
            ),
        )
    end

    CSV.write(TRANSITION_AUDIT_SUMMARY, DataFrame(summary_rows))
    CSV.write(TRANSITION_AUDIT_PER_TRANSITION, DataFrame(per_transition_rows))
    CSV.write(TRANSITION_AUDIT_WORST_SAMPLES, DataFrame(worst_rows))
    CSV.write(TRANSITION_AUDIT_GLOBAL_REASONS, DataFrame(global_rows))
    CSV.write(TRANSITION_AUDIT_SIMULATORS, DataFrame(simulator_rows))
    CSV.write(joinpath(TRANSITION_AUDIT_DIR, "indexing_diagnostics.csv"), DataFrame(indexing_rows_all))

    _audit_write_plots(per_transition_rows, global_rows, simulator_rows)
    _audit_write_diagnostics(summary_rows, indexing_rows_all)

    println("transition audit summary: ", TRANSITION_AUDIT_SUMMARY)
    println("transition audit per transition: ", TRANSITION_AUDIT_PER_TRANSITION)
    println("transition audit worst samples: ", TRANSITION_AUDIT_WORST_SAMPLES)
    println("global stress reasons: ", TRANSITION_AUDIT_GLOBAL_REASONS)
    println("simulator comparison: ", TRANSITION_AUDIT_SIMULATORS)
    println("plots: ", TRANSITION_AUDIT_PLOTS_DIR)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
