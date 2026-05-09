function build_double_pendulum_problem(cfg; pendulum_module)
    return pendulum_module.optimal_control_problem(;
        l1 = cfg.l1,
        l2 = cfg.l2,
        m1 = cfg.m1,
        m2 = cfg.m2,
        g = cfg.g,
        objective = cfg.objective,
    )
end

function build_double_pendulum_system_cfg(cfg; pendulum_module)
    problem = build_double_pendulum_problem(cfg; pendulum_module = pendulum_module)
    return (;
        concrete_system = problem.system,
        params = (; l1 = cfg.l1, l2 = cfg.l2, m1 = cfg.m1, m2 = cfg.m2, g = cfg.g),
    )
end

function build_double_pendulum_control_cfg(cfg; pendulum_module)
    problem = build_double_pendulum_problem(cfg; pendulum_module = pendulum_module)
    x0 = SVector(UT.get_center(problem.initial_set)...)
    return (;
        x0,
        initial_set = problem.initial_set,
        target_set = problem.target_set,
        problem,
    )
end

function build_double_pendulum_input_mapping(cfg)
    inputs = [[Float64(u)] for u in cfg.input_values]
    return MP.ListMapping(inputs)
end

function build_double_pendulum_symbolic_builder(cfg; pendulum_module)
    return function (prob, candidate, certifier_cfg)
        o = certifier_cfg.options
        return pendulum_module.symbolic_system(
            prob.system.X;
            l1 = cfg.l1,
            l2 = cfg.l2,
            m1 = cfg.m1,
            m2 = cfg.m2,
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
