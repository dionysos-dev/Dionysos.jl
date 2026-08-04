module TestMain

# Headless GR: render animations without a display (also needed under CI).
ENV["GKSwstype"] = "100"

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots   # loads the DionysosPlotsExt extension (provides animate_trajectory_dashboard)

# Smoke test for `animate_trajectory_dashboard` (DionysosPlotsExt): drive its 1D/2D state,
# 1D/2D input, hybrid-mode and error branches. Rendering is headless; we only assert that a
# non-nothing animation object (or the written file) comes back.
@testset "animate_trajectory_dashboard" begin
    xs = [SVector(0.0, 0.0), SVector(1.0, 0.5), SVector(2.0, 1.0)]
    us = [SVector(0.5, 0.2), SVector(0.4, 0.1)]
    traj = ST.Trajectory(xs; inputs = us)

    draw2d(fig, x, u) = scatter!(fig, [x[1]], [x[2]]; label = "")

    # 2D state / 2D input.
    anim =
        Dionysos.animate_trajectory_dashboard(draw2d, traj; xdims = (1, 2), udims = (1, 2))
    @test anim !== nothing

    # 1D state / 1D input (time-series panels).
    anim1 = Dionysos.animate_trajectory_dashboard(draw2d, traj; xdims = (1,), udims = (1,))
    @test anim1 !== nothing

    # frame_step > 1 renders a subset of frames.
    animf = Dionysos.animate_trajectory_dashboard(
        draw2d,
        traj;
        xdims = (1, 2),
        udims = (1, 2),
        frame_step = 2,
    )
    @test animf !== nothing

    # Hybrid trajectory: a `modes` channel adds the mode panel and passes the mode as a
    # 4th argument to `system_plot!`.
    htraj = ST.Trajectory(xs; inputs = us, modes = [1, 2, 2])
    draw_mode(fig, x, u, k) = scatter!(fig, [x[1]], [x[2]]; label = "")
    animm = Dionysos.animate_trajectory_dashboard(
        draw_mode,
        htraj;
        xdims = (1, 2),
        udims = (1, 2),
    )
    @test animm !== nothing

    # Writing to a .gif returns the filename.
    gifpath = joinpath(mktempdir(), "dash.gif")
    out = Dionysos.animate_trajectory_dashboard(
        draw2d,
        traj;
        xdims = (1, 2),
        udims = (1, 2),
        filename = gifpath,
    )
    @test out == gifpath
    @test isfile(gifpath)

    # Error branches: a trajectory without inputs, and an unsupported output extension.
    @test_throws ErrorException Dionysos.animate_trajectory_dashboard(
        draw2d,
        ST.Trajectory(xs);
        xdims = (1, 2),
    )
    @test_throws ErrorException Dionysos.animate_trajectory_dashboard(
        draw2d,
        traj;
        xdims = (1, 2),
        udims = (1, 2),
        filename = "dashboard.txt",
    )
end

end # module TestMain
