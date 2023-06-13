using JuMP
import HiGHS
mip_solver = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
@static if false
    # Does not pass test_Q_reuse
    import COSMO
    qp_solver = optimizer_with_attributes(
        COSMO.Optimizer,
        "max_iter" => 30000,
        MOI.Silent() => true,
    )
end
import OSQP
qp_solver = optimizer_with_attributes(
    OSQP.Optimizer,
    "polish" => 1,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true,
)
@static if false
    import Gurobi
    env = Gurobi.Env()
    qp_solver = optimizer_with_attributes(
        () -> Gurobi.Optimizer(env),
        MOI.Silent() => true,
        # Without this, I get NumericalError in BemporadMorari called by BranchAndBound
        "DualReductions" => 0,
    )
end
import Ipopt
cont_solver = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
import Pavito
miqp_solver = optimizer_with_attributes(
    Pavito.Optimizer,
    "mip_solver" => mip_solver,
    # Can use OSQP after https://github.com/jump-dev/Pavito.jl/pull/36
    "cont_solver" => cont_solver,
    MOI.Silent() => true,
)
