# The disturbance slide, run through the JuMP front-end exactly as shown:
# the same double integrator solved twice, once without `w` and once with it
# declared a DISTURBANCE, and the two winning sets drawn on top of each other.
# Output: presentation/robust_vs_nominal.png (overwrites)
#
# Run: julia --project=test presentation/asset_robust_jump.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP, Plots
import LazySets
using Symbolics, MathOptSymbolicAD
using Dionysos
const DI = Dionysos
const MP = DI.Mapping
const SY = DI.Symbolic

initial = LazySets.Hyperrectangle(; low = [-1.6, -0.1], high = [-1.4, 0.1])
target = LazySets.Hyperrectangle(; low = [-0.4, -0.4], high = [0.4, 0.4])
const w̄ = [0.02, 0.3]

function solve_it(; with_disturbance::Bool)
    model = Model(Dionysos.Optimizer)
    @variable(model, -2 <= x[1:2] <= 2)
    @variable(model, -1 <= u <= 1)
    if with_disturbance
        @variable(model, -w̄[i] <= w[i = 1:2] <= w̄[i])
        set_role!(w, Dionysos.DISTURBANCE)
        @constraint(model, ∂(x[1]) == x[2] + w[1])
        @constraint(model, ∂(x[2]) == u + w[2])
    else
        @constraint(model, ∂(x[1]) == x[2])
        @constraint(model, ∂(x[2]) == u)
    end
    @constraint(model, x in Start(initial))
    @constraint(model, x in Final(target))

    set_attribute(model, "time_step", 0.3)
    set_attribute(model, "jacobian_bound", u -> SMatrix{2, 2}(0.0, 0.0, 1.0, 0.0))
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "print_level", 0)
    optimize!(model)
    println("  termination: ", termination_status(model))

    win = get_attribute(model, "controllable_set")
    abs_sys = get_attribute(model, "abstract_system")
    return (win, SY.get_state_mapping(abs_sys))
end

println("nominal:")
nominal = solve_it(; with_disturbance = false)
println("robust:")
robust = solve_it(; with_disturbance = true)
println("nominal states: ", MP.get_n_state(nominal[1], nominal[2]))
println("robust  states: ", MP.get_n_state(robust[1], robust[2]))

fig = plot(; aspect_ratio = :equal, legend = :topright, dpi = 200)
plot!(
    fig,
    LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0]);
    color = :white,
    linecolor = :black,
    label = "",
)
plot!(
    fig,
    nominal;
    color = :steelblue,
    opacity = 0.5,
    linecolor = :steelblue,
    label = "nominal winning set",
)
plot!(
    fig,
    robust;
    color = :orange,
    opacity = 0.8,
    linecolor = :orange,
    label = "robust winning set (∀w)",
)
plot!(fig, target; fillalpha = 0.0, linecolor = :green, linewidth = 3, label = "target")
out = joinpath(@__DIR__, "robust_vs_nominal.png")
savefig(fig, out)
println("saved ", out)
