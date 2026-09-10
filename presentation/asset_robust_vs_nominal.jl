# Presentation asset: robust vs nominal winning sets — the same solver call, with
# and without a declared disturbance; the robust set is strictly inside the nominal.
# Output: presentation/robust_vs_nominal.png
#
# Run: julia --project=test presentation/asset_robust_vs_nominal.jl

ENV["GKSwstype"] = "100"

using StaticArrays
import LazySets
import MathematicalSystems
import MathOptInterface as MOI
using Plots
import Dionysos
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const AB = Dionysos.Optim.Abstraction

# A double integrator: position, velocity, bounded acceleration. The disturbance
# eats into the braking authority, so the robust winning set visibly loses the
# high-speed corners (the classic parabolic boundaries move inward).
X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
U = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
# ±0.3 on the acceleration channel: enough to carve a visible rim off the
# nominal winning set; ±0.4 already collapses the robust set onto the target
# (the folded inflation compounds along the braking arc).
W = LazySets.Hyperrectangle(; low = [-0.02, -0.3], high = [0.02, 0.3])

f_nom(x, u) = SVector(x[2], u[1])
f_noisy(x, u, w) = SVector(x[2] + w[1], u[1] + w[2])
jacobian_bound(u) = SMatrix{2, 2}(0.0, 0.0, 1.0, 0.0)

nominal = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(f_nom, 2, 1, X, U)
noisy = MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem(
    f_noisy,
    2,
    1,
    2,
    X,
    U,
    W,
)

initial_set = LazySets.Hyperrectangle(; low = [-1.6, -0.1], high = [-1.4, 0.1])
target_set = LazySets.Hyperrectangle(; low = [-0.4, -0.4], high = [0.4, 0.4])

function winning(system)
    problem = PR.OptimalControlProblem(
        system,
        initial_set,
        target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.25)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
    abs_sys = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    return (set, SY.get_state_mapping(abs_sys))
end

win_nom = winning(nominal)
win_rob = winning(noisy)
println("nominal winning states: ", MP.get_n_state(win_nom[1], win_nom[2]))
println("robust  winning states: ", MP.get_n_state(win_rob[1], win_rob[2]))

fig = plot(; aspect_ratio = :equal, legend = :topright, dpi = 200)
plot!(fig, X; color = :white, linecolor = :black, label = "")
plot!(
    fig,
    win_nom;
    color = :steelblue,
    opacity = 0.5,
    linecolor = :steelblue,
    label = "nominal winning set",
)
plot!(
    fig,
    win_rob;
    color = :orange,
    opacity = 0.8,
    linecolor = :orange,
    label = "robust winning set (∀w)",
)
plot!(fig, target_set; fillalpha = 0.0, linecolor = :green, linewidth = 3, label = "target")
out = joinpath(@__DIR__, "robust_vs_nominal.png")
savefig(fig, out)
println("saved ", out)
