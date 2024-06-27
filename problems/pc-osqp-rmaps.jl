module PCOSQPRM
using Test

using Dionysos: Dionysos
const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem

using FillArrays
using Polyhedra
using MathematicalSystems, HybridSystems
using SemialgebraicSets
using LinearAlgebra

using CDDLib: CDDLib

function system(lib, T::Type)

    # Generate the automaton nodes
    possible_values = [[-1, 0, 1], [-1, 0, 1], [-1, 0, 1], [0, 1], [0, 1], [0, 1]]
    nvars = length(possible_values)

    ## generate all possible combinations of the values
    function generate_combinations(possible_values, depth)
        if depth == 1
            return [[v] for v in possible_values[1]]
        end
        return [
            [v; c...] for v in possible_values[1] for c in generate_combinations(possible_values[2:end], depth - 1)
        ]
    end

    ## generate all possible combinations of the values
    nodes = generate_combinations(possible_values, nvars)
    nnodes = length(nodes)
    
	domains = []

	#Automaton:
	automaton = GraphAutomaton(nnodes)
    indexB = Int[]

	#Resetmaps: (guards + reset)
	A = T[
		0.9999981424449595 9.91220860939534e-12 3.4481569821196496e-7 9.205526977551878e-5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		-9.91220860939534e-12 0.9999981424449595 -9.205526977551878e-5 3.4481569821196496e-7 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		2.1535327064414276e-7 -2.644213073804883e-12 0.9999999080158245 -2.455696808386769e-5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		2.644213073804883e-12 2.1535327064414276e-7 2.455696808386769e-5 0.9999999080158245 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.9999999996930378 -2.477749999746475e-5 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 2.477749999746475e-5 0.9999999996930378 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.99875 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0012499999999999734 0.99875 0.0
		0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
	]

	B = T[
		6.256252307659102e-5 -3.128126153811649e-5 -3.128126153847453e-5 0.0 0.0 0.0
		-2.0671108273679513e-16 5.418073430928138e-5 -5.418073430907468e-5 0.0 0.0 0.0
		6.7365306173215874e-12 -3.368357184796113e-12 -3.3681734325413444e-12 0.0 0.0 0.0
		1.0608933096199221e-16 5.833953603352774e-12 -5.8340596926953345e-12 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0
		1.0 0.0 0.0 0.0 0.0 0.0
		0.0 1.0 0.0 0.0 0.0 0.0
		0.0 0.0 1.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.013888888888888593 0.013888888888888593 0.013888888888888593
		0.0 0.0 0.0 0.0 0.0 0.0
		0.0 0.0 0.0 0.0 0.0 0.0
	]

    C = T[
        1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
        0.0 1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -5.5 5.5 
    ]

	#ResetMaps = Vector{ConstrainedLinearControlDiscreteSystem}(undef, length(domains)^2)
	##cU = polyhedron(HalfSpace(T[0, 0, -1], 2) ∩ HalfSpace(T[0, 0, 1], 2), lib)

	function is_enabled(from, to)
        from_node = nodes[from]
        to_node = nodes[to]

        # On the first 3 elements of the nodes, we cannot go from -1 to 1 or 1 to -1
        for i in 1:3
            if from_node[i] == -1 && to_node[i] == 1
                return false
            end
            if from_node[i] == 1 && to_node[i] == -1
                return false
            end
        end

        return true
	end

	ntransitions = 0
	function maybe_add_transition(from, to)
        #global ntransitions
		if is_enabled(from, to)
			ntransitions += 1
			add_transition!(automaton, from, to, ntransitions)
            append!(indexB, from)
		end
	end
	for from in 1:nnodes, to in 1:nnodes
		maybe_add_transition(from, to)
	end

    modes = [ConstrainedContinuousIdentitySystem(12, i) for i in 1:nnodes]
    resetmaps = [AffineContinuousSystem(A, B*nodes[indexB[i]]) for i in 1:ntransitions]
	
	system = HybridSystem(
		automaton,
		# Modes
		modes,
		# Reset maps
		resetmaps,
        # Switching
		Fill(ControlledSwitching(), ntransitions),
	)

	system.ext[:q_T] = nmodes(system)
    system.ext[:q_C] = C
    
	return system
end

""""
	problem(lib, T::Type, gamma=0.95, x_0 = [1.0, -6.0], N = 11, tail_cost::Bool = false)

This function create the system with `PCOSQPRM.system`.

Then, we define initial conditions (continuous and discrete states) to this system
and set `N` as horizon, i.e., the number of allowed time steps.

We instantiate our Optimal Control Problem by defining the state and tail costs.
Notice that `state_cost` is defined to be \sum(\gamma^t \times norm(Cx_t)) for each time step `t`
and the `tail_cost` is defined to be \gamma^N V(x_N) that we considering constant in this version.

Notice that we used `Fill` for all `N` time steps as we consider time-invariant costs.
"""
function problem(
	lib = CDDLib.Library(),
	T::Type = Float64;
	gamma = 0.95,
	x_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	N = 3,
	tail_cost::Bool = false,
)
	sys = system(lib, T)
	if tail_cost
		transition_cost = Fill(UT.ZeroFunction(), nmodes(sys))
    end 

	state_cost = [
			mode == sys.ext[:q_T] ? UT.ConstantFunction(zero(T)) :
			UT.ConstantFunction(one(T)) for mode in modes(sys)
	]

	transition_cost = UT.QuadraticControlFunction(ones(T, 1, 1))
	problem = PR.OptimalControlProblem(
		sys,
		(x_0),
		sys.ext[:q_T],
		Fill(state_cost, N),
		Fill(Fill(transition_cost, ntransitions(sys)), N),
		N,
	)
	return problem
end

end

