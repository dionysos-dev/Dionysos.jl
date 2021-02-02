include("../gol_lazar_belta.jl")
include("../../test/solvers.jl")
using LinearAlgebra, Test
import CDDLib
using Dionysos

using Plots
using Colors
import GLPK
function text_in_set_plot!(pl, po, t, solver; kws...)
    plot!(pl, po; kws...)
    if t !== nothing
        c, r = hchebyshevcenter(hrep(po), solver, verbose=0)
        annotate!(pl, [(c..., text(t, 12))])
    end 
end    

    
function _prob( N, q0, x0, zero_cost::Bool)
    system = gol_lazar_belta(CDDLib.Library())
    if zero_cost
        state_cost = Fill(ZeroFunction(), nmodes(system))
    else
        state_cost = [mode == system.ext[:q_T] ? ConstantFunction(0.0) : ConstantFunction(1.0)
                        for mode in modes(system)]
    end
    return OptimalControlProblem(
        system,
        q0, x0,
        Fill(state_cost, N),
        Fill(Fill(QuadraticControlFunction(ones(1, 1)), ntransitions(system)), N),
        system.ext[:q_T],
        N
    )
end
function gol_lazar_belta_plot(system::HybridSystem, solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.ON))
    # Modes:
    p = Plots.plot(fmt = :png, fillcolor = :white)
    for mode in states(system)
        t = (system.ext[:q_T] in [mode, mode + 11]) ? "T" : (mode == system.ext[:q_A] ? "A" : (mode == system.ext[:q_B] ? "B" :
                mode <= 11 ? string(mode) : string(mode - 11)))
        text_in_set_plot!(p, stateset(system, mode), t, solver, fillcolor = :white, linecolor = :black)
    end

    # Obstacles:
    for i in eachindex(system.ext[:obstacles])
        text_in_set_plot!(p, system.ext[:obstacles][i], "O$i", solver, fillcolor = :black, fillalpha = 0.1)
    end

    return p
end



println("testing depth 9")
#_test(9, 8, [1.5, -2.5], false, true)



println("testing depth 11")
#_test(11, 3, [1.0, -6.0], false, true)
N = 11
q0 = 3
x0 = [1.0, -6.0]
zero_cost = false

# Pavito does not support indicator constraints yet so we use `false` here
algo = optimizer_with_attributes(BemporadMorari.Optimizer, 
"continuous_solver" => qp_solver, 
"mixed_integer_solver" => miqp_solver,
"indicator" => false, 
"log_level" => 0)

println("starting test")
problem = _prob(N, q0, x0, zero_cost)
println("Instantiating")
optimizer = MOI.instantiate(algo)
MOI.set(optimizer, MOI.RawParameter("problem"), problem)
println("Solving")
@time MOI.optimize!(optimizer)
println("Solved")
MOI.get(optimizer, MOI.TerminationStatus()) 
xu = MOI.get(optimizer, ContinuousTrajectoryAttribute())
#Initial state
#scatter!(p, [x0[1]], [x0[2]])
#annotate!(p, [(x0[1], x0[2] - 0.5, text("x0", 10))])

#Forward trajectories:

x1 = [xu.x[j][1] for j in eachindex(xu.x)]
x2 = [xu.x[j][2] for j in eachindex(xu.x)]

println(x1)
println(x2)
xu.u 
MOI.get(optimizer, MOI.ObjectiveValue()) 

if optimizer isa BranchAndBound.Optimizer
    optimizer.num_done + optimizer.num_pruned_bound + optimizer.num_pruned_inf 
    return optimizer.Q_function
    
end
p = gol_lazar_belta_plot(problem.system)
plot!(p, x1, x2)
display(p)
x=1
