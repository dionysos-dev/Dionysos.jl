module AStarAbstractionBuilder

using DataStructures: PriorityQueue
import Dionysos
const DO = Dionysos.Domain
const ST = Dionysos.System
const SY = Dionysos.Symbolic
const UT = Dionysos.Utils
const PR = Dionysos.Problem

"""
    build_abstraction_astar!(
        Xdom,
        AutomatonConstructor,
        concrete_problem;
        heuristic,
        synthesize_transition,
        synthesize_new_state,
        candidate_selector,
        concrete_system,
        max_iterations,
        unique_by_pos
    )

Entry point to construct a backward A*-based abstraction for a given problem.
"""
function build_abstraction_astar!(
    Xdom,
    AutomatonConstructor,
    concrete_problem;
    heuristic::Function,
    synthesize_transition::Function,
    synthesize_new_state::Function,
    candidate_selector::Function,
    concrete_system,
    max_iterations=10^5,
    unique_by_pos=true
)
    model = SY.LazyGridBasedSymbolicModel(Xdom, AutomatonConstructor)

    setI = concrete_problem.initial_set[1]
    setT = concrete_problem.target_set[1]

    posI = DO.get_pos_by_coord(Xdom, UT.get_center(setI))
    posT = DO.get_pos_by_coord(Xdom, UT.get_center(setT))

    qI = add_state!(model, posI, setI)
    qT = add_state!(model, posT, setT)

    abstract_problem = PR.OptimalControlProblem(model, [qI], [qT], nothing, nothing)

    return build_abstraction_astar!(
        model,
        qI,
        qT,
        heuristic,
        synthesize_transition,
        synthesize_new_state,
        candidate_selector,
        concrete_system,
        abstract_problem;
        max_iterations,
        unique_by_pos
    )
end

"""
    build_abstraction_astar!(model, qI, qT, ...)

Core A*-like backward abstraction builder.
"""
function build_abstraction_astar!(
    model::SY.LazyGridBasedSymbolicModel,
    qI::Int,
    qT::Int,
    heuristic::Function,
    synthesize_transition::Function,
    synthesize_new_state::Function,
    candidate_selector::Function,
    concrete_system,
    abstract_problem;
    max_iterations=10^5,
    unique_by_pos=true
)
    frontier = PriorityQueue{Int, Float64}()
    visited = Set{Int}()
    grid = DO.get_grid(SY.get_state_domain(model))
    Xdom = SY.get_state_domain(model)
    initial_pos = DO.get_pos_by_state(model, qI)

    push!(frontier, qT => 0.0)  # root of the abstraction tree

    iteration = 0
    while !isempty(frontier) && iteration < max_iterations
        iteration += 1

        q_curr = dequeue!(frontier)
        pos_curr = get_pos(model, q_curr)
        cost_curr = SY.get_best_worst_case_cost(model, q_curr, abstract_problem)

        # 1. Select nearby positions
        candidate_positions = candidate_selector(grid, pos_curr)
        candidate_positions = filter(pos -> pos ∈ Xdom, candidate_positions)

        for pos in candidate_positions
            center = DO.get_coord_by_pos(Xdom, pos)
            h_val = heuristic(center)

            q = nothing
            κ = nothing
            cost = Inf

            # 2. Try reusing an existing abstract state
            reuse_possible = false
            if unique_by_pos || pos == initial_pos
                existing_states = get_states_by_pos(model, pos)
                if !isempty(existing_states)
                    q = existing_states[1]
                    κ, cost = synthesize_transition(center, get_concrete_states(model, q_curr), concrete_system)
                    reuse_possible = true
                end
            end

            # 3. If not reusable, create new abstract state
            if !reuse_possible || κ === nothing
                set, κ, cost = synthesize_new_state(center, get_concrete_states(model, q_curr), concrete_system)
                if set === nothing || κ === nothing
                    continue  # skip invalid transitions
                end
                q = add_state!(model, pos, set)
            end

            s = add_input!(model, κ)
            SY.add_transition!(model, q, q_curr, s)
            SY.add_transition_cost!(model, q, s, cost)

            # 4. A* cost propagation
            g_val = cost_curr + cost
            if q ∉ visited
                push!(visited, q)
                frontier[q] = g_val + h_val
            else
                prev_cost = SY.get_best_worst_case_cost(model, q, abstract_problem)
                if g_val < prev_cost
                    frontier[q] = g_val + h_val
                end
            end

            # 5. Early termination if initial reached
            if q == qI
                println("A* abstraction reached initial set at pos = $pos, state = $q")
                return
            end
        end
    end

    println("A* abstraction terminated after $iteration iterations.")
end

end  # module AStarAbstractionBuilder
