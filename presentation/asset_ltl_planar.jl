# Single-panel animation of a co-safe LTL run: the labelled regions as a fixed
# background, the trajectory unrolling over them. The full dashboard wastes two
# thirds of the frame on an integrator, whose "system view" carries nothing.
# Output: presentation/gallery/{integrator_ltl,integrator_ltl_until}_planar.gif
#
# Run: julia --project=test presentation/asset_ltl_planar.jl [seq|until]

ENV["GKSwstype"] = "100"

which = isempty(ARGS) ? "seq" : ARGS[1]
root = dirname(@__DIR__)

if which == "seq"
    include(joinpath(root, "examples", "Integrator", "integrator_ltl_1.jl"))
    colors = Dict(:g1 => :red, :g2 => :cyan, :g3 => :orange, :obs => :black)
    out = joinpath(@__DIR__, "gallery", "integrator_ltl_planar.gif")
else
    include(joinpath(root, "examples", "Integrator", "integrator_ltl_2.jl"))
    colors = Dict(:g1 => :red, :g2 => :cyan, :danger => :orange, :obs => :black)
    out = joinpath(@__DIR__, "gallery", "integrator_ltl_until_planar.gif")
end

using Plots
xs = [Float64(s[1]) for s in ST.states(traj)]
ys = [Float64(s[2]) for s in ST.states(traj)]

anim_planar = Plots.Animation()
ks = unique(vcat(collect(2:2:length(xs)), length(xs)))
for k in ks
    f = plot(;
        aspect_ratio = :equal,
        legend = false,
        size = (760, 700),
        dpi = 150,
        title = string(φ),
        titlefontsize = 9,
    )
    plot!(f, concrete_problem; ap_colors = colors, label = false)
    plot!(f, xs[1:k], ys[1:k]; color = :blue, linewidth = 2, label = false)
    scatter!(f, [xs[k]], [ys[k]]; color = :blue, markersize = 5, label = false)
    Plots.frame(anim_planar, f)
end
## hold on the finished picture
for _ in 1:10
    Plots.frame(anim_planar, current())
end

gif(anim_planar, out; fps = 8)
println("saved ", out)
