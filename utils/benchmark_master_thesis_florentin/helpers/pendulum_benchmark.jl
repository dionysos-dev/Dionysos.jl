function build_pendulum_problem(cfg; pendulum_module)
    return pendulum_module.optimal_control_problem(;
        l = cfg.l,
        g = cfg.g,
        objective = cfg.objective,
    )
end

function build_pendulum_system_cfg(cfg; pendulum_module)
    problem = build_pendulum_problem(cfg; pendulum_module = pendulum_module)
    return (;
        concrete_system = problem.system,
        params = (; l = cfg.l, g = cfg.g),
    )
end

function build_pendulum_control_cfg(cfg; pendulum_module)
    problem = build_pendulum_problem(cfg; pendulum_module = pendulum_module)
    x0 = SVector(UT.get_center(problem.initial_set)...)
    return (;
        x0,
        initial_set = problem.initial_set,
        target_set = problem.target_set,
        problem,
    )
end

function build_pendulum_input_mapping(cfg)
    inputs = [[Float64(u)] for u in cfg.input_values]
    return MP.ListMapping(inputs)
end

function build_pendulum_symbolic_builder(cfg; pendulum_module)
    return function (prob, candidate, certifier_cfg)
        o = certifier_cfg.options
        return pendulum_module.symbolic_system(
            prob.system.X;
            l = cfg.l,
            g = cfg.g,
            _U_ = prob.system.U,
            Ts = candidate.Ts,
            ΔX = o.ΔX,
            ΔU = o.ΔU,
            ΔW = o.ΔW,
            rk4_num_substeps = o.symbolic_rk4_substeps,
        )
    end
end

function build_pendulum_wrap_state(cfg)
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

function periodic_state_error(x, xref, cfg)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    if x isa SVector{N} where {N}
        return typeof(x)(Tuple(e))
    end
    return e
end

pendulum_state_error(x, xref, cfg) = periodic_state_error(x, xref, cfg)

function project_input_to_domain(u, u_domain)
    uvec = collect(Float64, u)
    u_domain isa UT.HyperRectangle || return uvec
    return collect(clamp.(uvec, u_domain.lb, u_domain.ub))
end

project_pendulum_input_to_domain(u, u_domain) = project_input_to_domain(u, u_domain)
