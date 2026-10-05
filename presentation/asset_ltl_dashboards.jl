# Re-animate the two co-safe LTL drivers with the labelled regions drawn in the
# state-space panel (the dashboard's `state_background!` hook), so the formula
# can be read off the picture.
# Output: presentation/gallery/{integrator,pendulum}_ltl.gif
#
# Run: julia --project=test presentation/asset_ltl_dashboards.jl <which>
#   <which> = "integrator" (default) or "pendulum"

ENV["GKSwstype"] = "100"

which = isempty(ARGS) ? "integrator" : ARGS[1]
root = dirname(@__DIR__)

if which == "integrator"
    include(joinpath(root, "examples", "Integrator", "integrator_ltl_1.jl"))
    colors = Dict(:g1 => :red, :g2 => :cyan, :g3 => :orange, :obs => :black)
    out = joinpath(@__DIR__, "gallery", "integrator_ltl.gif")
    dims = (1, 2)
    xlab, ylab = "x₁", "x₂"
elseif which == "until"
    include(joinpath(root, "examples", "Integrator", "integrator_ltl_2.jl"))
    colors = Dict(:g1 => :red, :g2 => :cyan, :danger => :orange, :obs => :black)
    out = joinpath(@__DIR__, "gallery", "integrator_ltl_until.gif")
    dims = (1, 2)
    xlab, ylab = "x₁", "x₂"
    system_plot! = Integrator.system_plot!()
else
    include(joinpath(root, "examples", "Pendulum", "simple_pendulum_co_safe_ltl.jl"))
    colors = Dict(:g1 => :red, :g2 => :cyan, :obs => :black)
    out = joinpath(@__DIR__, "gallery", "pendulum_ltl.gif")
    dims = (1, 2)
    xlab, ylab = "θ [rad]", "ω [rad/s]"
end

using Plots

## the labelled regions, drawn under the trajectory at every frame
function regions_background!(p_state, x)
    plot!(p_state, concrete_problem; ap_colors = colors, label = false)
    return p_state
end

anim_regions = Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = dims,
    udims = (1,),
    frame_step = 2,
    xlabel_state = xlab,
    ylabel_state = ylab,
    state_background! = regions_background!,
)
gif(anim_regions, out; fps = 6)
println("saved ", out)
