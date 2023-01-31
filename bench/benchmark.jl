PROBLEMS_PATH = joinpath(@__DIR__, "..", "problems")
all_problems = []
all_solvers = []
for file in readdir(PROBLEMS_PATH)
    if file == "robot_problem.jl"
        continue
    end
    path = joinpath(PROBLEMS_PATH, file)
    if endswith(file, ".jl")
        mod = include(path)
        if !Base.isdefined(mod, :problem)
            @warn("$mod defined in $path/$file does not define a problem function, skipping this benchmark.")
        end
        @info("Loading problem from $mod in $path/$file")
        @time push!(all_problems, mod.problem())
    end
end

import OSQP
osqp = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true
);

import Ipopt
ipopt = optimizer_with_attributes(
    Ipopt.Optimizer,
    MOI.Silent() => true
);

QP_SOLVERS = [("OSQP", osqp), ("Ipopt", ipopt)]

import HiGHS
highs = optimizer_with_attributes(
    HiGHS.Optimizer,
    MOI.Silent() => true
);

MIP_SOLVERS = [("HiGHS", highs)]

MIQP_SOLVERS = []
import Pavito
for (mip_name, mip_solver) in MIP_SOLVERS
    for (qp_name, qp_solver) in QP_SOLVERS
        pavito = optimizer_with_attributes(
            Pavito.Optimizer,
            "mip_solver" => mip_solver,
            "cont_solver" => qp_solver,
            MOI.Silent() => true
        )
        push!(MIQP_SOLVERS, ("Pavito($mip_name,$qp_name)", pavito))
    end
end

SOLVERS = []

for (qp_name, qp_solver) in QP_SOLVERS
    for (miqp_name, miqp_solver) in MIQP_SOLVERS
        bemporad_morari = optimizer_with_attributes(BemporadMorari.Optimizer{Float64},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0
        )
        push!(all_solvers, ("BemporadMorary($qp_name,$miqp_name)", bemporad_morari))
    end
end

function bench(solver, problem)
    optimizer = MOI.instantiate(solver)
    try
        MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
    catch err
        if err isa MOI.UnsupportedAttribute
            return
        else
            rethrow(err)
        end
    end
    MOI.optimize!(optimizer)

    # We check the solver time
    time = MOI.get(optimizer, MOI.SolveTimeSec())
    @show time

    # the termination status
    termination = MOI.get(optimizer, MOI.TerminationStatus())
    @show termination
    return
end

for (name, solver) in all_solvers
    println(name)
    for problem in all_problems
        bench(solver, problem)
    end
end
