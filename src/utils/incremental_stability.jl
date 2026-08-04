# Grid sizing for center-simulation abstractions of incrementally stable (δ-GAS) systems.
#
# When a uniform-grid abstraction only simulates cell *centers*, the abstraction is unsound
# on its own: nothing bounds how far trajectories starting elsewhere in a cell drift from
# the simulated center. For an incrementally stable (δ-GAS) system a common quadratic
# Lyapunov function bounds that drift, and the bound both sizes the grid and gives the
# precision ε of the resulting approximately-bisimilar symbolic model
# [girard2010approximately](@cite).
#
# For a linear switched system the δ-GAS certificate *is* the common quadratic Lyapunov
# function `V(x, y) = ‖x - y‖_P`: the error `e = x - y` obeys the same switched dynamics, so
# the stability certificate is also the incremental one. Obtain `P` and the per-step
# contraction `γ` (the factor by which `V` shrinks each sampling step) from
# `compute_quadratic_pieces_pclf` applied to the sampled modes `exp(A_σ · τ)` with a
# single-node graph; the helpers below turn `(P, γ)` into a grid step or a precision.

"""
    grid_step_from_lyapunov(P, ε, γ; state_dim = size(P, 1)) -> η

Largest uniform grid step `η` for which a **center-simulation** uniform-grid abstraction of
a δ-GAS system is `ε`-approximately bisimilar to the sampled system, following
[girard2010approximately](@cite).

`P` is the symmetric positive-definite matrix of a common quadratic Lyapunov function
`V(x, y) = ‖x - y‖_P`, `γ ∈ [0, 1)` is its per-step contraction (`V` shrinks by a factor `γ`
each sampling step — for a continuous-time rate `κ` sampled at `τ`, `γ = e^{-κτ}`), and `ε`
is the target precision in the Euclidean (output) norm.

The abstraction relation `{(x, x_q) : ‖x - x_q‖_P ≤ √λ_min(P) · ε}` is preserved when the
per-step quantization error is absorbed by one step of contraction:

    √λ_max(P) · (η √n / 2) ≤ (1 - γ) · √λ_min(P) · ε,

so

    η = 2 (1 - γ) · √(λ_min(P) / λ_max(P)) · ε / √n .

The condition number of `P` therefore directly controls how fine the grid must be.

Obtain `(P, γ)` from a common quadratic Lyapunov function, e.g. `compute_quadratic_pieces_pclf`
on the sampled modes. See [`precision_from_grid_step`](@ref) for the inverse.
"""
function grid_step_from_lyapunov(P, ε, γ; state_dim::Int = size(P, 1))
    0 <= γ < 1 ||
        error("per-step contraction γ must lie in [0, 1); the system must be contractive.")
    Ps = LinearAlgebra.Symmetric(Matrix{Float64}(P))
    λmin = LinearAlgebra.eigmin(Ps)
    λmax = LinearAlgebra.eigmax(Ps)
    return 2 * (1 - γ) * ε * sqrt(λmin / λmax) / sqrt(state_dim)
end

"""
    precision_from_grid_step(P, η, γ; state_dim = size(P, 1)) -> ε

Precision `ε` (Euclidean norm) of the `ε`-approximately bisimilar center-simulation
abstraction obtained with uniform grid step `η`; the inverse of
[`grid_step_from_lyapunov`](@ref):

    ε = η √n · √(λ_max(P) / λ_min(P)) / (2 (1 - γ)) .
"""
function precision_from_grid_step(P, η, γ; state_dim::Int = size(P, 1))
    0 <= γ < 1 ||
        error("per-step contraction γ must lie in [0, 1); the system must be contractive.")
    Ps = LinearAlgebra.Symmetric(Matrix{Float64}(P))
    λmin = LinearAlgebra.eigmin(Ps)
    λmax = LinearAlgebra.eigmax(Ps)
    return η * sqrt(state_dim) * sqrt(λmax / λmin) / (2 * (1 - γ))
end
