# The pipeline slide's left panel: one figure carrying the three stages of the
# pipeline as three trajectories in the (θ, ω) plane.
#
#   1. the grid-abstraction seed        (Plan)
#   2. the MPPI refinement of that seed (Refine)
#   3. the closed loop under the certified funnel controller, started off the
#      nominal trajectory on the boundary of the entry funnel (Certify)
#
# It runs the real demo and reads its results, so the picture and the numbers on
# the slide come from the same run.
#
# Output: presentation/gallery/pendulum_three_traj.png (overwrites)
#
# Run: julia --project=bench presentation/asset_pendulum_three_traj.jl

ENV["GKSwstype"] = "100"

include(
    joinpath(
        dirname(@__DIR__),
        "research",
        "TrajectoryCertificationOptimizer",
        "demo_pendulum.jl",
    ),
)

import LinearAlgebra

# The plant map the certificate was built on: one RK4 step of size Δt.
function rk4_step(x, u)
    g = SimplePendulum.dynamic(params)
    k1 = collect(g(x, u))
    k2 = collect(g(x .+ (Δt / 2) .* k1, u))
    k3 = collect(g(x .+ (Δt / 2) .* k2, u))
    k4 = collect(g(x .+ Δt .* k3, u))
    return SVector{2}(x .+ (Δt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4))
end

# A start on the boundary of the entry funnel: the certificate says every such
# point is recovered, so this is the trajectory the certificate is about.
function funnel_start(ctrl)
    E = ctrl.ellipsoids[1]
    Q = Matrix(LazySets.shape_matrix(E))
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Q))
    axis = F.vectors[:, argmax(F.values)] .* sqrt(maximum(F.values))
    return SVector{2}(collect(LazySets.center(E)) .+ 0.98 .* axis)
end

function closed_loop(ctrl, x0)
    xs = [x0]
    x = x0
    for k in 1:length(ctrl.kappas)
        u = ST.output_control(ctrl, k, x)
        u === nothing && break
        x = rk4_step(x, u)
        push!(xs, x)
    end
    return xs
end

mppi_xs = collect(ST.states(lifted_x))

# `lift` shifts each trajectory onto the branch where its OWN endpoint sits near
# the target, which can land the seed a full period away from the refinement.
# Put it on the refinement's branch instead, so the two are comparable.
function on_branch_of(traj, reference_end)
    xs = collect(ST.states(ST.unwrap_trajectory(traj, (1,), (2π,))))
    shift = 2π * round((xs[end][1] - reference_end[1]) / (2π))
    return [SVector(x[1] - shift, x[2]) for x in xs]
end

seed_xs = on_branch_of(seed_traj, mppi_xs[end])

let raw = collect(ST.states(seed_traj))[end]
    println("seed endpoint (wrapped)   : ", round.(raw; digits = 3))
    println("seed endpoint (on branch) : ", round.(seed_xs[end]; digits = 3))
    println("mppi endpoint             : ", round.(mppi_xs[end]; digits = 3))
    println("seed in target            : ", wrap(raw) ∈ problem.target_set)
end

fig = plot(;
    xlabel = "θ  [rad]",
    ylabel = "ω  [rad/s]",
    title = "Certified pendulum swing-up",
    legend = :topleft,
    size = (760, 620),
    dpi = 160,
)

plot!(fig, problem.initial_set; color = :gray, alpha = 0.5, label = "initial set")
plot!(fig, problem.target_set; color = :green, alpha = 0.35, label = "target set")

for (i, E) in enumerate(funnel)
    plot!(
        fig,
        E;
        color = :steelblue,
        alpha = 0.22,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = i == 1 ? "certified funnel" : "",
    )
end

plot!(
    fig,
    [x[1] for x in seed_xs],
    [x[2] for x in seed_xs];
    color = :gray40,
    linewidth = 2,
    linestyle = :dash,
    label = "1. grid abstraction seed",
)

plot!(
    fig,
    [x[1] for x in mppi_xs],
    [x[2] for x in mppi_xs];
    color = :black,
    linewidth = 2.5,
    label = "2. MPPI refinement",
)

# Drawn last and dotted: it tracks the refinement closely, which is the point,
# so a solid line on top would simply hide it.
if loop_success
    x0 = funnel_start(ctrl)
    cl_xs = closed_loop(ctrl, x0)
    println("closed loop: $(length(cl_xs)) steps from $(round.(x0; digits = 3))")
    plot!(
        fig,
        [x[1] for x in cl_xs],
        [x[2] for x in cl_xs];
        color = :darkorange,
        linewidth = 3,
        linestyle = :dot,
        label = "3. closed loop, off-nominal start",
    )
end

out = joinpath(@__DIR__, "gallery", "pendulum_three_traj.png")
savefig(fig, out)
println("saved ", out)
