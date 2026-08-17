# Lyapunov state feedback for affine systems. Unrelated to the transition
# synthesis SDPs — no ellipsoids, no noise, no cost bound — just the classical
# discrete-time stabilizability LMI.

import LazySets
using JuMP

"""
    stabilizing_feedback(subsys, sdp_solver) -> (feasible, K, P, γ)

For a stabilizable affine system, find the state-feedback gain `K` and the
matrix `P` satisfying the discrete-time Lyapunov inequality
`(A+BK)'P(A+BK) − P ≺ 0`, minimizing the condition number of `P` (`γ` is the
attained smallest eigenvalue of `P⁻¹`). On an unsuccessful solve the result is
`(false, nothing, nothing, nothing)` — no values are extracted from a failed
model.
"""
function stabilizing_feedback(
    subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem,
    sdp_solver,
)
    A = subsys.A
    B = subsys.B
    n = size(A, 1)
    m = size(B, 2)

    model = Model(sdp_solver)
    @variable(model, L[i = 1:m, j = 1:n])
    @variable(model, S[i = 1:n, j = 1:n], PSD)
    @variable(model, gamma)

    id = Matrix{Float64}(LA.I, n, n)

    @constraint(
        model,
        [
            S transpose(A * S + B * L)
            A * S+B * L S
        ] >= 1e-10 * Matrix{Float64}(LA.I, 2n, 2n),
        PSDCone()
    )
    @constraint(model, id >= S, PSDCone())
    @constraint(model, S >= gamma * id, PSDCone())

    @objective(model, Max, gamma)

    optimize!(model)

    termination_status(model) == MOI.OPTIMAL || return false, nothing, nothing, nothing

    P = inv(value.(S))
    K = value.(L) * P
    return true, K, P, value(gamma)
end
