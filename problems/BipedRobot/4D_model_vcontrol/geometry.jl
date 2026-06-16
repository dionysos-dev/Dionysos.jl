module Geometry

using StaticArrays

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

"""
Get the cartesian coordinates of each joint from the four angular coordinates and the foot on the ground
"""
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


"""
Get the 4 angular coordinates from the cartesian coordinates of the joints
"""
function get_angular_coordinates(cartesian_coords::Tuple, geometry)

    lh, lk, lf, rh, rk, rf = cartesian_coords
    grounded_left_foot = (cartesian_coords[3] == SVector(0.0, 0.0))

    theta1 = asin((lk[1]-lh[1])/geometry.hip_to_knee) - geometry.theta_hip
    theta2 = asin((lf[1]-lk[1])/geometry.knee_to_foot) - geometry.theta_hip - theta1
    theta3 = asin((rk[1]-rh[1])/geometry.hip_to_knee) - geometry.theta_hip
    theta4 = asin((rf[1]-rk[1])/geometry.knee_to_foot) - geometry.theta_hip - theta3

    return SVector(theta1, theta2, theta3, theta4), grounded_left_foot
end

"""
Useful for visualization
"""
function robot_segments(theta, geometry, grounded_left_foot)

    lh, lk, lf, rh, rk, rf =
        get_cartesian_coordinates(theta, geometry, grounded_left_foot)

    left_x = [lf[1], lk[1], lh[1]]
    left_y = [lf[2], lk[2], lh[2]]

    right_x = [rf[1], rk[1], rh[1]]
    right_y = [rf[2], rk[2], rh[2]]

    return left_x, left_y, right_x, right_y
end

end