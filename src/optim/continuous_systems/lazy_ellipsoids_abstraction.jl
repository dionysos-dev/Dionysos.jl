export LazyEllipsoidsAbstraction

module LazyEllipsoidsAbstraction
using JuMP
import LinearAlgebra as LA

import MathematicalSystems
MS = MathematicalSystems

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

"""
    Optimizer{T} <: Dionysos.Optim.AbstractDionysosOptimizer

Abstraction-based solver using a lazy abstraction method with ellipsoidal cells (RRT-based).
"""
mutable struct Optimizer{T} <: OP.AbstractDionysosOptimizer
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, UT.Tree}
    abstract_system_full::Union{Nothing, Any}
    abstract_controller::Union{Nothing, Any}
    concrete_controller::Union{Nothing, ST.AbstractContinuousController}
    abstract_lyap_fun::Union{Nothing, Any}
    concrete_lyap_fun::Union{Nothing, Any}

    # algorithm parameters
    distance::Union{Nothing, Any}
    rand_state::Union{Nothing, Any}
    new_conf::Union{Nothing, Any}
    keep::Union{Nothing, Any}
    stop_crit::Union{Nothing, Any}
    RRTstar::Union{Nothing, Bool}
    compute_transition::Union{Nothing, Any}
    maxIter::Union{Nothing, Int}
    maxδx::Union{Nothing, Any}
    maxδu::Union{Nothing, Any}
    λ::Union{Nothing, Any}
    sdp_opt::Union{Nothing, Any}
    k1::Union{Nothing, Int}
    k2::Union{Nothing, Int}
    continues::Union{Nothing, Bool}

    _found_initial::Bool
    _can_rewire::Bool
    _initial_node::Union{Nothing, Any}

    solve_time_sec::T

    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            false,
            nothing,
            0.0,
        )
    end
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "concrete_problem" && !(value isa PR.OptimalControlProblem)
        throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
    end
    return OP.set_field_attribute!(model, param, value)
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function set_optimizer!(
    opt::Optimizer,
    concrete_problem,
    sdp_opt,
    maxδx,
    maxδu,
    λ,
    k1,
    k2,
    RRTstar,
    continues,
    maxIter,
)
    # reset internal state
    opt._found_initial = false
    opt._can_rewire = false
    opt._initial_node = nothing

    # user params
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("sdp_opt"), sdp_opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxδx"), maxδx)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxδu"), maxδu)
    MOI.set(opt, MOI.RawOptimizerAttribute("λ"), λ)
    MOI.set(opt, MOI.RawOptimizerAttribute("k1"), k1)
    MOI.set(opt, MOI.RawOptimizerAttribute("k2"), k2)
    MOI.set(opt, MOI.RawOptimizerAttribute("RRTstar"), RRTstar)
    MOI.set(opt, MOI.RawOptimizerAttribute("continues"), continues)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxIter"), maxIter)

    # defaults (set here so your old call sites still work)
    MOI.set(opt, MOI.RawOptimizerAttribute("distance"), distance)
    MOI.set(opt, MOI.RawOptimizerAttribute("rand_state"), rand_state)
    MOI.set(opt, MOI.RawOptimizerAttribute("new_conf"), new_conf)
    MOI.set(opt, MOI.RawOptimizerAttribute("keep"), keep)
    MOI.set(opt, MOI.RawOptimizerAttribute("stop_crit"), stop_crit)
    MOI.set(opt, MOI.RawOptimizerAttribute("compute_transition"), compute_transition)

    return
end

function set_optimizer!(
    opt::Optimizer,
    concrete_problem,
    sdp_opt,
    maxδx,
    maxδu,
    λ,
    k1,
    k2,
    RRTstar,
    continues,
    maxIter,
    distance,
    rand_state,
    new_conf,
    keep,
    stop_crit,
    compute_transition,
)
    opt._found_initial = false
    opt._can_rewire = false
    opt._initial_node = nothing

    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("sdp_opt"), sdp_opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxδx"), maxδx)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxδu"), maxδu)
    MOI.set(opt, MOI.RawOptimizerAttribute("λ"), λ)
    MOI.set(opt, MOI.RawOptimizerAttribute("k1"), k1)
    MOI.set(opt, MOI.RawOptimizerAttribute("k2"), k2)
    MOI.set(opt, MOI.RawOptimizerAttribute("RRTstar"), RRTstar)
    MOI.set(opt, MOI.RawOptimizerAttribute("continues"), continues)
    MOI.set(opt, MOI.RawOptimizerAttribute("maxIter"), maxIter)

    MOI.set(opt, MOI.RawOptimizerAttribute("distance"), distance)
    MOI.set(opt, MOI.RawOptimizerAttribute("rand_state"), rand_state)
    MOI.set(opt, MOI.RawOptimizerAttribute("new_conf"), new_conf)
    MOI.set(opt, MOI.RawOptimizerAttribute("keep"), keep)
    MOI.set(opt, MOI.RawOptimizerAttribute("stop_crit"), stop_crit)
    MOI.set(opt, MOI.RawOptimizerAttribute("compute_transition"), compute_transition)

    return
end

# ----------------------------
# Controller / Lyapunov helpers
# ----------------------------

struct TreeStaticController{TR} <: ST.AbstractContinuousController
    tree::TR
end

ST.controller_kind(::TreeStaticController) = ST.StaticKind()
ST.domain(ctrl::TreeStaticController) = ctrl.tree
ST.initial_state(::TreeStaticController) = nothing
ST.update_state(::TreeStaticController, q, x) = nothing

function _best_tree_action(ctrl::TreeStaticController, x)
    compare(E, x) = x ∈ E
    nodes = UT.get_nodes(ctrl.tree, x, compare)
    isempty(nodes) && return nothing

    sorted_nodes = sort(nodes; by = UT.compare)
    isempty(sorted_nodes) && return nothing

    local_map = UT.get_action(sorted_nodes[1])
    local_map === nothing && return nothing

    u = MS.apply(local_map, x)
    u === nothing && return nothing

    return u
end

function ST.is_defined(ctrl::TreeStaticController, q, x)
    return _best_tree_action(ctrl, x) !== nothing
end

function ST.output_control(ctrl::TreeStaticController, q, x)
    return _best_tree_action(ctrl, x)
end

function build_concrete_controller(abstract_tree::UT.Tree)
    return TreeStaticController(abstract_tree)
end

function build_abstract_lyap_fun()
    return abstract_lyap_fun(node) = UT.get_path_cost(node)
end

function build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)
    compare(E, x) = x ∈ E
    function concrete_lyap_fun(x)
        nodes = UT.get_nodes(abstract_system, x, compare)
        sorted_nodes = sort(nodes; by = UT.compare)
        isempty(sorted_nodes) && return Inf
        return abstract_lyap_fun(sorted_nodes[1])
    end
    return concrete_lyap_fun
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    concrete_problem = optimizer.concrete_problem
    optimizer.abstract_problem = concrete_problem

    # Co-design abstract system and abstract controller
    abstract_system = UT.RRT(
        concrete_problem.target_set,
        concrete_problem.initial_set,
        optimizer.distance,
        optimizer.rand_state,
        optimizer.new_conf,
        optimizer.keep,
        optimizer.stop_crit,
        optimizer;
        maxIter = optimizer.maxIter,
        RRTstar = optimizer.RRTstar,
        compute_transition = optimizer.compute_transition,
        k1 = optimizer.k1,
        k2 = optimizer.k2,
    )

    optimizer.abstract_system = abstract_system

    # Build concrete controller and Lyapunov function
    optimizer.concrete_controller = build_concrete_controller(abstract_system)
    abstract_lyap_fun = build_abstract_lyap_fun()
    optimizer.abstract_lyap_fun = abstract_lyap_fun
    optimizer.concrete_lyap_fun =
        build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)

    optimizer.solve_time_sec = time() - t_ref
    return
end

# ----------------------------
# Default algorithm parameters
# ----------------------------

function distance(E1::UT.Ellipsoid, E2::UT.Ellipsoid)
    return UT.centerDistance(E1, E2)
end

function get_candidate(
    tree::UT.Tree,
    X,
    E0::UT.Ellipsoid;
    probTarget = 0.15,
    probSkew = 0.35,
)
    guess = UT.sample(X)
    r = rand()

    if r < probTarget
        return E0.c
    elseif r < probTarget + probSkew
        α = 0.7 + 0.3 * rand()   # strongly biased toward E0
        return α * E0.c + (1 - α) * guess
    else
        return guess
    end
end

function rand_state(
    tree::UT.Tree,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    opt::Optimizer,
)
    concrete_problem = opt.concrete_problem
    xrand = get_candidate(tree, concrete_problem.system.X, EI)
    return UT.Ellipsoid(Matrix{Float64}(LA.I(length(xrand))), xrand)
end

function get_closest_reachable_point(
    concrete_system,
    xinit,
    xtarget,
    U,
    Uformat;
    nSamples = 500,
)
    unew = UT.sample(U)
    wnew = zeros(concrete_system.nw)
    xnew = concrete_system.f_backward_eval(xinit, unew, wnew)
    uBestDist = LA.norm(xnew - xtarget)

    for i in 1:nSamples
        ucandnew = UT.sample(U) * 0.002 * i
        xcandnew = concrete_system.f_backward_eval(xinit, ucandnew, wnew)
        if LA.norm(xcandnew - xtarget) < uBestDist
            uBestDist = LA.norm(xcandnew - xtarget)
            xnew = xcandnew
            unew = ucandnew
        end
    end

    return (unew, xnew, uBestDist)
end

function new_conf(
    abstract_system::UT.Tree,
    Nnear::UT.NodeT,
    Erand::UT.Ellipsoid,
    opt::Optimizer,
)
    concrete_problem = opt.concrete_problem
    concrete_system = concrete_problem.system

    (unew, xnew, uBestDist) = get_closest_reachable_point(
        concrete_system,
        Nnear.state.c,
        Erand.c,
        concrete_system.U,
        concrete_system.Uformat,
    )

    wnew = zeros(concrete_system.nw)

    (affineSys, L) = ST.buildAffineApproximation(
        concrete_system.fsymbolic,
        concrete_system.x,
        concrete_system.u,
        concrete_system.w,
        xnew,
        unew,
        wnew,
        xnew .+ concrete_system.ΔX,
        unew .+ concrete_system.ΔU,
        wnew .+ concrete_system.ΔW,
    )

    S = UT.get_full_psd_matrix(concrete_problem.transition_cost)

    return UT.transition_backward(
        affineSys,
        Nnear.state,
        xnew,
        unew,
        concrete_system.Uformat,
        concrete_system.Wformat,
        S,
        L,
        opt.sdp_opt;
        λ = opt.λ,
        maxδx = opt.maxδx,
        maxδu = opt.maxδu,
    )
end

function keep(
    abstract_system::UT.Tree,
    LSACnew,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    opt::Optimizer;
    scale_for_obstacle = true,
)
    concrete_problem = opt.concrete_problem
    obstacles = concrete_problem.system.obstacles

    minDist = Inf
    iMin = 0

    for (i, data) in enumerate(LSACnew)
        Enew, cont, cost, Nnear = data

        if Enew === nothing
            # infeasible candidate
        elseif EI ∈ Enew
            iMin = i
            break
        elseif minDist > LA.norm(EI.c - Enew.c)
            if Nnear == abstract_system.root || LA.eigmin(EI.P * 0.5 - Enew.P) > 0
                iMin = i
                minDist = LA.norm(EI.c - Enew.c)
            end
        end
    end

    iMin == 0 && return []

    ElMin, contMin, costMin, NnearMin = LSACnew[iMin]

    if ElMin !== nothing
        if all(O -> !UT.is_intersecting(ElMin, O), obstacles)
            return [LSACnew[iMin]]
        elseif scale_for_obstacle
            for O in obstacles
                ElMin = UT.compress_if_intersection(ElMin, O)
                ElMin === nothing && return []
            end
            return [(ElMin, contMin, costMin, NnearMin)]
        else
            return []
        end
    end

    return []
end

function compute_transition(E1::UT.Ellipsoid, E2::UT.Ellipsoid, opt::Optimizer)
    concrete_problem = opt.concrete_problem
    concrete_system = concrete_problem.system

    xnew = E1.c
    unew = zeros(concrete_system.nu)
    wnew = zeros(concrete_system.nw)

    (affineSys, L) = ST.buildAffineApproximation(
        concrete_system.fsymbolic,
        concrete_system.x,
        concrete_system.u,
        concrete_system.w,
        xnew,
        unew,
        wnew,
        xnew .+ concrete_system.ΔX,
        unew .+ concrete_system.ΔU,
        wnew .+ concrete_system.ΔW,
    )

    S = UT.get_full_psd_matrix(concrete_problem.transition_cost)

    ans, cont, cost = UT.transition_fixed(
        affineSys,
        E1,
        E2,
        concrete_system.Uformat,
        concrete_system.Wformat,
        S,
        opt.sdp_opt,
    )

    return ans, cont, cost
end

function stop_crit(
    abstract_system::UT.Tree,
    LNnew,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    opt::Optimizer,
)
    continues = opt.continues
    minDist = 10.0

    for Nnew in LNnew
        E = Nnew.state
        if distance(EI, E) <= minDist
            ans, cont, cost = compute_transition(EI, E, opt)
            if ans
                if !opt._found_initial
                    opt._initial_node = UT.add_node!(abstract_system, EI, Nnew, cont, cost)
                    println("Path cost from EI : ", opt._initial_node.path_cost)
                    opt._found_initial = true
                end

                if !continues
                    return true
                else
                    if opt._can_rewire &&
                       (cost + Nnew.path_cost < opt._initial_node.path_cost)
                        UT.rewire(abstract_system, opt._initial_node, Nnew, cont, cost)
                        println("Path cost from EI : ", opt._initial_node.path_cost)
                    end
                end
                opt._can_rewire = true
            end
        end
    end

    newEllipsoids = [newNode.state for newNode in LNnew]
    return any(E -> (EI ∈ E), newEllipsoids) && !continues
end

end # module
