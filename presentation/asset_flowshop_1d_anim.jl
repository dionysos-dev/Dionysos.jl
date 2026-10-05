# Animated version of the FlowShop 1D state-over-time figure: the guard windows
# and final target stay fixed, the trajectory unrolls, switches appear as they
# happen, and the animation holds on the complete picture.
# Output: presentation/gallery/flowshop_1d_anim.gif
#
# Run: julia --project=test presentation/asset_flowshop_1d_anim.jl

ENV["GKSwstype"] = "100"

include(
    joinpath(
        dirname(@__DIR__),
        "examples",
        "FlowShopScheduling",
        "flow_shop_scheduling_1D.jl",
    ),
)

# Rebuild the driver's background (guards + per-task grids + final target),
# reusing the variables the include left in Main.
function background()
    bg = plot(;
        title = "1D flowshop scheduling",
        xlabel = "Time (s)",
        ylabel = "System state (x)",
        xlims = (t_min - t_margin, t_max + t_margin),
        ylims = (x_min - x_margin, x_max + x_margin),
        legend = :outerright,
        legendtitle = "Legend",
        size = (1400, 900),
        dpi = 150,
    )
    region_colors = palette(:tab20)
    for (i, t) in enumerate(HS.transitions(concrete_system.automaton))
        guard = HS.guard(concrete_system, t)
        guard === nothing && continue
        xa = [LazySets.low(guard, 1), LazySets.high(guard, 1)]
        ta = [LazySets.low(guard, 2), LazySets.high(guard, 2)]
        c = region_colors[(i - 1) % length(region_colors) + 1]
        plot!(
            bg,
            [ta[1], ta[2], ta[2], ta[1], ta[1]],
            [xa[1], xa[1], xa[2], xa[2], xa[1]];
            fill = (0, 0.15),
            fillcolor = c,
            alpha = 0.3,
            linecolor = c,
            linewidth = 2,
            linestyle = :dash,
            label = "Guard $(i): x∈[$(xa[1]),$(xa[2])], t∈[$(ta[1]),$(ta[2])]",
        )
    end
    xa = [LazySets.low(final_x_target, 1), LazySets.high(final_x_target, 1)]
    ta = [final_spec.tmin, final_spec.tmax]
    plot!(
        bg,
        [ta[1], ta[2], ta[2], ta[1], ta[1]],
        [xa[1], xa[1], xa[2], xa[2], xa[1]];
        fill = (0, 0.25),
        fillcolor = :magenta,
        alpha = 0.4,
        linecolor = :magenta,
        linewidth = 3,
        label = "Final target",
    )
    return bg
end

bg0 = background()

function overlay(k)
    f = deepcopy(bg0)
    plot!(
        f,
        t_traj[1:k],
        x_traj[1:k];
        color = :black,
        linewidth = 2,
        alpha = 0.8,
        label = "Trajectory",
    )
    for i in switch_indices
        i <= k || continue
        plot!(
            f,
            [t_traj[i - 1], t_traj[i]],
            [x_traj[i - 1], x_traj[i]];
            color = :red,
            linewidth = 3,
            linestyle = :dash,
            alpha = 0.9,
            label = false,
        )
    end
    sw = [i for i in switch_indices if i <= k]
    if !isempty(sw)
        scatter!(
            f,
            [t_traj[i] for i in sw],
            [x_traj[i] for i in sw];
            color = :red,
            marker = :diamond,
            markersize = 8,
            label = "Switch",
        )
    end
    scatter!(
        f,
        [t_traj[1]],
        [x_traj[1]];
        color = :green,
        marker = :star5,
        markersize = 8,
        label = "Initial state",
    )
    scatter!(f, [t_traj[k]], [x_traj[k]]; color = :black, markersize = 5, label = false)
    return f
end

anim = Plots.Animation()
N = length(t_traj)
ks = collect(2:2:N)
last(ks) == N || push!(ks, N)
for k in ks
    Plots.frame(anim, overlay(k))
end
final_frame = overlay(N)
for _ in 1:12
    Plots.frame(anim, final_frame)
end

out = joinpath(@__DIR__, "gallery", "flowshop_1d_anim.gif")
Plots.gif(anim, out; fps = 8)
println("saved ", out)
