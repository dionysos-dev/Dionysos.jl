using Plots

"""
To save the animation in gif format.
"""
function animate_robot(
    X,
    dt,
    geometry;
    grounded_left_foot=true,
    filename="./problems/BipedRobot/4D_model_vcontrol/walking_robot.gif"
)

    anim = @animate for x in X

        lx, ly, rx, ry =robot_segments(x, geometry, grounded_left_foot)

        p = plot(xlim=(-0.5, 0.5), ylim=(-0.05, 0.45), aspect_ratio=:equal, legend=false, xlabel="x (m)", ylabel="y (m)")

        # left leg
        plot!(p, lx, ly, lw=4, marker=:circle)
        # right leg
        plot!(p, rx, ry, lw=4, marker=:circle)
        # hip connection
        plot!(p, [lx[end], rx[end]], [ly[end], ry[end]], lw=4)

        p
    end

    gif(anim, filename, fps=div(100, Int(100*dt)))
end

"""
To visualize the animation in a live plot. It does not save the animation.
"""
function animate_robot_live(
    X,
    dt,
    geometry;
    grounded_left_foot=true
)

    for x in X
        lx, ly, rx, ry =robot_segments(x, geometry, grounded_left_foot)
        # left leg
        p = plot(lx, ly, lw=4, marker=:circle, aspect_ratio=:equal, xlim=(-0.5, 0.5), ylim=(-0.05, 0.45), legend=false, xlabel="x (m)", ylabel="y (m)")
        # right leg
        plot!(p, rx, ry, lw=4, marker=:circle)
        # hip connection
        plot!(p, [lx[end], rx[end]], [ly[end], ry[end]], lw=4)

        display(p)

        sleep(dt)
    end
end