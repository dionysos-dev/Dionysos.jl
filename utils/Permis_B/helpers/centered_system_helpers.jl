import StaticArrays: SVector

"""
    build_uniform_grid_abstraction(concrete_system, Δt, hx, Udom, jacobian_bound; ...)

Construit l'abstraction uniforme Dionysos (mode `CENTER_SIMULATION` par defaut).

- `hx`: pas de grille par dimension d'etat.
- `Udom`: alphabet d'actions discret.
- `jacobian_bound`: borne sur la croissance locale utilisee par l'abstraction.
- options periodiques: utiles pour les dimensions angulaires (repliement modulo `2π`).
"""
function build_uniform_grid_abstraction(
    concrete_system,
    Δt,
    hx,
    Udom,
    jacobian_bound;
    approx_mode = AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    # Configuration du solveur d'abstraction via attributs MOI.
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("Udom"), Udom)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), approx_mode)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

    if with_period
        MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periodic_periods)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
    end

    MOI.optimize!(optimizer)
    return optimizer
end

"""
    build_uniform_grid_controller!(optimizer, concrete_system, _I_, _T_)

Attache un probleme de controle optimal a l'abstraction puis synthétise
un controleur concret atteignant `_T_` depuis `_I_`.
"""
function build_uniform_grid_controller!(optimizer, concrete_system, _I_, _T_)
    concrete_problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)

    MOI.optimize!(optimizer)

    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    target_set = concrete_problem.target_set
    return concrete_controller, target_set
end

"""
    simulate_closed_loop(concrete_system, concrete_controller, Δt, x0, target_set; ...)

Simule la boucle fermee sur le systeme discretise.
Si `with_period=true`, les dimensions periodiques sont repliees avant:
- le test d'arret (`x ∈ target_set`),
- le stockage de la trajectoire.
"""
function simulate_closed_loop(
    concrete_system,
    concrete_controller,
    Δt,
    x0,
    target_set;
    nstep = 700,
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    disc = ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)

    periodic_wrapper =
        with_period ?
        ST.get_periodic_wrapper(periodic_dims, periodic_periods; start = periodic_start) :
        (x -> x)

    reached(x) = (periodic_wrapper(x) ∈ target_set)

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        disc,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
        wrap = periodic_wrapper,
    )
    return x_traj, u_traj
end

"""
    get_Udom(_U_, hu)

Construit un domaine d'entrées discret (`DomainList`) a partir d'un pas `hu`,
puis intersecte avec l'ensemble continu `_U_`.
"""
function get_Udom(_U_, hu)
    u0 = zeros(SVector{2, Float64})
    input_grid = DO.GridFree(u0, hu)
    Udom = Dionysos.Domain.DomainList(input_grid)
    Dionysos.Domain.add_set!(Udom, _U_, Dionysos.Domain.CENTER)
    return Udom
end

"""
    collect_trajectory_lists(x_traj, u_traj)

Convertit les trajectoires Dionysos (ou iterables generiques) en vecteurs Julia.
Pratique pour post-traitement, affichage et LMIs.
"""
function collect_trajectory_lists(x_traj, u_traj)
    state_list = x_traj isa ST.Trajectory ? copy(x_traj.seq) : collect(x_traj)
    input_list = u_traj isa ST.Trajectory ? copy(u_traj.seq) : collect(u_traj)
    return state_list, input_list
end
