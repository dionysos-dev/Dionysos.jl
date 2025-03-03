export LazyEllipsoidsAbstraction

module LazyEllipsoidsAbstraction
using LinearAlgebra, JuMP, IntervalArithmetic, Random

Random.seed!(0)

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

# temp
global myBool = true
global myBool2 = false
global NI = nothing

"""
    Optimizer{T} <: MOI.AbstractOptimizer

Abstraction-based solver using the lazy abstraction method with ellipsoidal cells.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, UT.Tree}
    abstract_system_full::Union{Nothing, Any}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Union{Nothing, Any}
    abstract_lyap_fun::Union{Nothing, Any}
    concrete_lyap_fun::Union{Nothing, Any}

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
            0.0,
        )
    end
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "concrete_problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function set_optimizer!(
    optimizer::Optimizer,
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
    global myBool = true
    global myBool2 = false
    global NI = nothing
    # User's paramaters
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_opt"), sdp_opt)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxδx"), maxδx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxδu"), maxδu)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("λ"), λ)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("k1"), k1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("k2"), k2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("RRTstar"), RRTstar)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("continues"), continues)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxIter"), maxIter)
    # Defaults's algorihthm pararaters
    MOI.set(optimizer, MOI.RawOptimizerAttribute("distance"), distance)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("rand_state"), rand_state)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("new_conf"), new_conf)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("keep"), keep)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("stop_crit"), stop_crit)
    return MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("compute_transition"),
        compute_transition,
    )
end

function set_optimizer!(
    optimizer::Optimizer,
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
    # User's paramaters
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_opt"), sdp_opt)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxδx"), maxδx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxδu"), maxδu)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("λ"), λ)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("k1"), k1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("k2"), k2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("RRTstar"), RRTstar)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("continues"), continues)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxIter"), maxIter)
    # Defaults's algorihthm pararaters
    MOI.set(optimizer, MOI.RawOptimizerAttribute("distance"), distance)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("rand_state"), rand_state)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("new_conf"), new_conf)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("keep"), keep)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("stop_crit"), stop_crit)
    return MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("compute_transition"),
        compute_transition,
    )
end

function build_concrete_controller(abstract_system)
    compare(E, x) = x ∈ E
    function concrete_controller(x)
        nodes = UT.get_nodes(abstract_system, x, compare)
        sorted_nodes = sort(nodes; by = UT.compare)
        isempty(sorted_nodes) && return nothing
        cont = UT.get_action(sorted_nodes[1])
        c_eval = ST.get_c_eval(cont)
        return c_eval(x)
    end
    return concrete_controller
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
    # Co-design the abstract system and the abstract controller
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
    # Build the concrete controller
    optimizer.concrete_controller = build_concrete_controller(abstract_system)
    abstract_lyap_fun = build_abstract_lyap_fun()
    optimizer.abstract_lyap_fun = abstract_lyap_fun
    optimizer.concrete_lyap_fun =
        build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)

    optimizer.solve_time_sec = time() - t_ref
    return
end

# ### Optimizer's parameters

# SI, SF
function distance(E1::UT.Ellipsoid, E2::UT.Ellipsoid)
    return UT.centerDistance(E1, E2)
end

function get_candidate(
    tree::UT.Tree,
    X::IntervalBox,
    E0::UT.Ellipsoid;
    probSkew = 0.0,
    probE0 = 0.05,
    intialDist = 1,
)
    guess = UT.sample(X)
    randVal = rand()
    if randVal > probSkew + probE0
        return guess
    elseif randVal > probSkew
        return E0.c
    else
        closestNode, dist =
            UT.findNClosestNode(tree, UT.Ellipsoid(Matrix{Float64}(I(nx)), E0))
        l = randVal / probSkew
        r = dist / intialDist
        return (E0.c * l + closestNode.state.c * (1 - l)) * (1 - 0.3 * r) +
               (0.3 * r) * guess
    end
end

function rand_state(
    tree::UT.Tree,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    optimizer::Optimizer,
)
    concrete_problem = optimizer.concrete_problem
    xrand = get_candidate(tree, concrete_problem.system.X, EI)
    return UT.Ellipsoid(Matrix{Float64}(I(length(xrand))), xrand)
end

# data-driven technique on nominal system (without noise)
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
    uBestDist = norm(xnew - xtarget)
    for i in 1:nSamples
        ucandnew = UT.sample(U) * 0.002 * i
        xcandnew = concrete_system.f_backward_eval(xinit, ucandnew, wnew)
        if norm(xcandnew - xtarget) < uBestDist
            uBestDist = norm(xcandnew - xtarget)
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
    optimizer::Optimizer,
)
    concrete_problem = optimizer.concrete_problem
    concrete_system = concrete_problem.system
    (unew, xnew, uBestDist) = get_closest_reachable_point(
        concrete_system,
        Nnear.state.c,
        Erand.c,
        concrete_system.U,
        concrete_system.Uformat,
    )
    wnew = zeros(concrete_system.nw)
    X̄ = IntervalBox(xnew .+ concrete_system.ΔX)
    Ū = IntervalBox(unew .+ concrete_system.ΔU)
    W̄ = IntervalBox(wnew .+ concrete_system.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        concrete_system.fsymbolic,
        concrete_system.x,
        concrete_system.u,
        concrete_system.w,
        xnew,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )
    S = UT.get_full_psd_matrix(concrete_problem.transition_cost)
    return SY.transition_backward(
        affineSys,
        Nnear.state,
        xnew,
        unew,
        concrete_system.Uformat,
        concrete_system.Wformat,
        S,
        L,
        optimizer.sdp_opt;
        λ = optimizer.λ,
        maxδx = optimizer.maxδx,
        maxδu = optimizer.maxδu,
    )
end

# heuristic: keep only the closest ellipsoid to the initial ellipsoid
function keep(
    abstract_system::UT.Tree,
    LSACnew,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    optimizer::Optimizer;
    scale_for_obstacle = true,
)
    concrete_problem = optimizer.concrete_problem
    obstacles = concrete_problem.system.obstacles
    minDist = Inf
    iMin = 0
    for (i, data) in enumerate(LSACnew)
        Enew, cont, cost, Nnear = data
        if Enew === nothing
            print("\tInfeasible")
        elseif EI ∈ Enew
            iMin = i
            break
        elseif minDist > norm(EI.c - Enew.c) # minPathCost > cost + Eclosest.path_cost
            if Nnear == abstract_system.root || eigmin(EI.P * 0.5 - Enew.P) > 0 # E ⊂ E0 => P-P0>0
                iMin = i
                minDist = norm(EI.c - Enew.c)
            else
            end
        else
        end
    end
    if iMin == 0
        return []
    end
    ElMin, contMin, costMin, NnearMin = LSACnew[iMin]
    if ElMin !== nothing
        if all(O -> !UT.is_intersected(ElMin, O), obstacles)
            return [LSACnew[iMin]]
        elseif scale_for_obstacle
            for O in obstacles
                ElMin = UT.compress_if_intersection(ElMin, O)
                if ElMin === nothing
                    return []
                end
            end
            return [(ElMin, contMin, costMin, NnearMin)]
        else
            return []
        end
    else
        return []
    end
end

function compute_transition(E1::UT.Ellipsoid, E2::UT.Ellipsoid, optimizer::Optimizer)
    concrete_problem = optimizer.concrete_problem
    concrete_system = concrete_problem.system
    xnew = E1.c
    unew = zeros(concrete_system.nu)
    wnew = zeros(concrete_system.nw)
    X̄ = IntervalBox(xnew .+ concrete_system.ΔX)
    Ū = IntervalBox(unew .+ concrete_system.ΔU)
    W̄ = IntervalBox(wnew .+ concrete_system.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        concrete_system.fsymbolic,
        concrete_system.x,
        concrete_system.u,
        concrete_system.w,
        xnew,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )
    S = UT.get_full_psd_matrix(concrete_problem.transition_cost)
    ans, cont, cost = SY.transition_fixed(
        affineSys,
        E1,
        E2,
        concrete_system.Uformat,
        concrete_system.Wformat,
        S,
        optimizer.sdp_opt,
    )
    return ans, cont, cost
end

function stop_crit(
    abstract_system::UT.Tree,
    LNnew,
    EF::UT.Ellipsoid,
    EI::UT.Ellipsoid,
    distance,
    optimizer::Optimizer,
)
    continues = optimizer.continues
    minDist = 10.0
    for Nnew in LNnew
        E = Nnew.state
        if distance(EI, E) <= minDist
            ans, cont, cost = compute_transition(EI, E, optimizer)
            if ans
                if myBool
                    global NI = UT.add_node!(abstract_system, EI, Nnew, cont, cost)
                    println("Path cost from EI : ", NI.path_cost)
                    global myBool = false
                end
                if !continues
                    return true
                else
                    if myBool2 && cost + Nnew.path_cost < NI.path_cost
                        UT.rewire(abstract_system, NI, Nnew, cont, cost)
                        println("Path cost from EI : ", NI.path_cost)
                    end
                end
                global myBool2 = true
            end
        end
    end
    newEllipsoids = [newNode.state for newNode in LNnew]
    return any(map(E -> (EI ∈ E), newEllipsoids)) && !continues
end

end
