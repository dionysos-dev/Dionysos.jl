module PostProc

using StaticArrays
using Plots
using LazySets
using ..Geometry


const PLOT = plot(aspect_ratio=:equal, xlim=(-0.5, 0.5), ylim=(-0.05, 0.45), legend=false, xlabel="x (m)", ylabel="y (m)")


"""
General function to visualize a robot configuration
"""

function plot_robot(
    x,
    geometry,
    grounded_left_foot
)

    lx, ly, rx, ry = Geometry.robot_segments(x, geometry, grounded_left_foot)
    # left leg
    p = plot(lx, ly, lw=4, marker=:circle, label="Left leg", aspect_ratio=:equal, xlim=(-0.5, 0.5), ylim=(-0.05, 0.45), legend=false, xlabel="x (m)", ylabel="y (m)")
    # right leg
    plot!(p, rx, ry, lw=4, marker=:circle, label="Right leg")
    # hip connection
    plot!(p, [lx[end], rx[end]], [ly[end], ry[end]], lw=4, label="")

    return p

end

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

        p = plot_robot(x, geometry, grounded_left_foot)

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
    grounded_left_foot=true,
    p = PLOT
)

    for x in X
        p = plot_robot(x, geometry, grounded_left_foot)
        display(p)
        sleep(dt)
    end

    return p
end

"""
Plot the target set in the 2D plane (propagates angle intervals, starting from the grounded foot)
"""
function plot_configuration_interval(
    theta_lb::SVector{N},
    theta_ub::SVector{N},
    geometry::Geometry.geometry,
    grounded_left_foot;
    p = PLOT
) where {N}

    for theta1 in theta_lb[1]:1e-1:theta_ub[1]
        for theta2 in theta_lb[2]:1e-1:theta_ub[2]
            for theta3 in theta_lb[3]:1e-1:theta_ub[3]
                for theta4 in theta_lb[4]:1e-1:theta_ub[4]
                    for point in Geometry.get_cartesian_coordinates(SVector(theta1, theta2, theta3, theta4), geometry, grounded_left_foot)
                        scatter!(p, [point[1]], [point[2]], markercolor=:blue, markeralpha=0.3, label="")
                    end
                end
            end
        end
    end

    display(p)

    return p

end

end