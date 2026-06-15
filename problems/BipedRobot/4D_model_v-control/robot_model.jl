using StaticArrays
using Plots

include("./postproc.jl")

"""
States (angles) and inputs (angular velocities):
    1 -> left hip
    2 -> left knee
    3 -> right hip
    4 -> right knee
"""

robot_dynamics(x,u,dt) = x + dt*u

"""
Robot geometry:
    Torso angle: 0°
    Hip-to-knee length: 20.2cm
    Knee-to-ankle length: 17.2cm
"""
struct geometry
    theta_hip
    hip_to_knee
    knee_to_foot
end

function get_cartesian_coordinates(theta::SVector, geometry, grounded_left_foot)
    (theta1, theta2, theta3, theta4) = grounded_left_foot ? (theta[1], theta[2], theta[3], theta[4]) : (theta[3], theta[4], theta[1], theta[2])

    grounded_foot = SVector(0.0, 0.0)
    grounded_knee = grounded_foot + geometry.knee_to_foot*SVector(-sin(geometry.theta_hip+theta1+theta2), cos(geometry.theta_hip+theta1+theta2))
    grounded_hip = grounded_knee + geometry.hip_to_knee*SVector(-sin(geometry.theta_hip+theta1), cos(geometry.theta_hip+theta1))
    swing_hip = grounded_hip
    swing_knee = swing_hip + geometry.hip_to_knee*SVector(sin(geometry.theta_hip+theta3), -cos(geometry.theta_hip+theta3))
    swing_foot = swing_knee + geometry.knee_to_foot*SVector(sin(geometry.theta_hip+theta3+theta4), -cos(geometry.theta_hip+theta3+theta4))

    return grounded_left_foot ? (grounded_hip, grounded_knee, grounded_foot, swing_hip, swing_knee, swing_foot) : (swing_hip, swing_knee, swing_foot, grounded_hip, grounded_knee, grounded_foot)
end

function simulate_robot(x0, U, dt)
    N = length(U)

    X = Vector{typeof(x0)}(undef, N + 1)
    X[1] = x0

    for k in 1:N
        X[k+1] = robot_dynamics(X[k], U[k], dt)
    end

    return X
end

function robot_segments(theta, geometry, grounded_left_foot)

    lh, lk, lf, rh, rk, rf =
        get_cartesian_coordinates(theta, geometry, grounded_left_foot)

    left_x = [lf[1], lk[1], lh[1]]
    left_y = [lf[2], lk[2], lh[2]]

    right_x = [rf[1], rk[1], rh[1]]
    right_y = [rf[2], rk[2], rh[2]]

    return left_x, left_y, right_x, right_y
end



"""
USAGE EXAMPLE
"""

robot_geometry = geometry(0, 0.202, 0.172)
x0 = SVector(0.0, 0.0, pi/4, -pi/4)

dt = 0.01
T = 5

lx, ly, rx, ry = robot_segments(x0, robot_geometry, true)

plot(
    lx, ly,
    lw=4,
    marker=:circle,
    aspect_ratio=:equal,
    xlim=(-0.5,0.5),
    ylim=(-0.05,0.45),
    legend=false
)

U = [
    SVector(
        0.0,
        0.0,
        1*sin(0.05*k),
        -1*sin(0.05*k)
    )
    for k in 1:N
]

X = simulate_robot(x0, U, dt)
animate_robot_live(
    X,
    dt,
    robot_geometry,
    grounded_left_foot=true
)