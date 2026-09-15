# Getting-started pendulum dashboard with the forbidden band drawn in both
# panels: as a wedge from the pivot in the pendulum view (left) and as a shaded
# vertical band in the state-space panel (top right), with the target in green.
# Output: presentation/gallery/pendulum_dashboard.gif (overwrites)
#
# Run: julia --project=test presentation/asset_pendulum_getting_started.jl

ENV["GKSwstype"] = "100"

include(joinpath(dirname(@__DIR__), "docs", "src", "examples", "getting_started.jl"))

band_lo = -π + 16π / 180
band_hi = -π + 38π / 180

base_plot! = SimplePendulum.system_plot!()
function pendulum_with_band!(fig, x, u)
    plot!(fig; aspect_ratio = :equal)
    θs = range(band_lo, band_hi; length = 24)
    xs = [0.0; [1.25 * sin(θ) for θ in θs]; 0.0]
    ys = [0.0; [-1.25 * cos(θ) for θ in θs]; 0.0]
    plot!(
        fig,
        xs,
        ys;
        seriestype = :shape,
        color = :black,
        opacity = 0.35,
        linecolor = :black,
        label = false,
    )
    return base_plot!(fig, x, u)
end

function band_background!(p_state, x)
    plot!(
        p_state,
        (winning_set, SY.get_state_mapping(abstract_system));
        color = :purple,
        opacity = 0.15,
        linecolor = :purple,
        linealpha = 0.15,
        label = false,
    )
    plot!(
        p_state,
        LazySets.Hyperrectangle(; low = [band_lo, -10.0], high = [band_hi, 10.0]);
        color = :black,
        opacity = 0.25,
        linecolor = :black,
        label = false,
    )
    plot!(
        p_state,
        UT.set_in_period(upright, SVector(1), SVector(2π), SVector(-π));
        color = :green,
        opacity = 0.3,
        linecolor = :green,
        label = false,
    )
    return p_state
end

anim2 = Dionysos.animate_trajectory_dashboard(
    pendulum_with_band!,
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.1,
    frame_step = 2,
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    ylabel_input = "τ [N·m]",
    state_background! = band_background!,
)

out = joinpath(@__DIR__, "gallery", "pendulum_dashboard.gif")
gif(anim2, out; fps = 6)
println("saved ", out)
