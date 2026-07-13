module RobotModel

using StaticArrays
# using Plots
using Dionysos
using MathematicalSystems
using LinearAlgebra

using ..Geometry
# using ..PostProc

const UT = Dionysos.Utils
const PB = Dionysos.Problem
const ST = Dionysos.System
const MP = Dionysos.Mapping

"""
States (angles) and inputs (angular velocities):
    1 -> left hip
    2 -> left knee
    3 -> right hip
    4 -> right knee
"""

robot_dynamics() = (x,u) -> SVector{4}(u[1], u[2], u[3], u[4])
robot_dynamics_dt(dt) = (x,u) -> x + dt*u

# For ẋ = u, ∂f/∂x = 0, so a constant zero matrix bound is valid.
# Returned function signature matches Dionysos' growth-bound constructors: L(u).
function jacobian_bound()
    return (u) -> @SMatrix [
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0
    ]
end

function system(;
    # theta_max = 30°
    # omega_max = 15°/s
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        robot_dynamics(),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function dt_system(
    dt;
    # theta_max = 30°
    # omega_max = 15°/s
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6))
)
    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        robot_dynamics_dt(dt),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(;
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6)),
    _I_ = UT.HyperRectangle(SVector(-pi/12, -pi/12, -pi/12, -pi/12), SVector(pi/12, pi/12, pi/12, pi/12)),
    _T_ = UT.HyperRectangle(SVector(-pi/3+pi/12, -pi/12, pi/3-pi/12, -pi/12), SVector(-pi/3, pi/12, pi/3, pi/12)),
)
    sys = system(; _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end



function index_vectors(lb::NTuple{N,Int}, ub::NTuple{N,Int}) where {N}
    ranges = ntuple(i -> lb[i]:ub[i], N)
    return (SVector{N}(t) for t in Iterators.product(ranges...))
end

function remove_infeasible_cells(_X_::UT.HyperRectangle{N, T} , state_grid::MP.GridFree, infeasible_shape::Any, geometry, grounded_left_foot) where {N, T}

    infeasible_sets = typeof(_X_)[]

    lb = MP.get_pos_by_coord(state_grid, _X_.lb)
    ub = MP.get_pos_by_coord(state_grid, _X_.ub)
    for idx in index_vectors(lb, ub)
        center = MP.get_coord_by_pos(state_grid, idx)
        rec = MP.get_elem_by_coord(state_grid, center)
        shrinked_rec = UT.HyperRectangle(
            rec.lb * 0.95,
            rec.ub * 0.95
        )
        cartesian = Geometry.get_cartesian_coordinates(center, geometry, grounded_left_foot)
        swing_foot_coords = grounded_left_foot ? cartesian[6] : cartesian[3]
        if Base.in(swing_foot_coords, infeasible_shape) || swing_foot_coords[2] < -1e-8
            push!(infeasible_sets, shrinked_rec)
        end
    end

    return UT.LazySetMinus(_X_, UT.LazySetUnion(infeasible_sets))

end


function compute_target_set(
    state_grid::MP.GridFree,
    swing_foot_coords::SVector{2},
    geometry,
    grounded_left_foot
)

    cell = MP.get_elem_by_coord(state_grid, SVector(0.0, 0.0, 0.0, 0.0))
    target_sets = typeof(cell)[]
    count = 0

    Lmax = geometry.hip_to_knee + geometry.knee_to_foot
    Lmin = abs(geometry.hip_to_knee - geometry.knee_to_foot)
    for hip_x in 0:1e-3:swing_foot_coords[1]
        for hip_y in swing_foot_coords[2]:1e-3:Lmax
            ds = norm([hip_x, hip_y])
            dw = norm([hip_x, hip_y] - swing_foot_coords)
            if Lmin <= ds <= Lmax && Lmin <= dw <= Lmax
                theta_k_stance = -acos(1 + (ds^2 - Lmax^2)/(2*geometry.hip_to_knee*geometry.knee_to_foot))

                theta_k_swing = -acos(1 + (dw^2 - Lmax^2)/(2*geometry.hip_to_knee*geometry.knee_to_foot))

                v = inv([-geometry.knee_to_foot*cos(theta_k_stance)-geometry.hip_to_knee -geometry.knee_to_foot*sin(theta_k_stance); -geometry.knee_to_foot*sin(theta_k_stance) geometry.knee_to_foot*cos(theta_k_stance)+geometry.hip_to_knee]) * [hip_x; hip_y]

                theta_h_stance = atan(v[1], v[2]) - geometry.theta_hip

                v = inv([-geometry.knee_to_foot*cos(theta_k_swing)-geometry.hip_to_knee -geometry.knee_to_foot*sin(theta_k_swing); -geometry.knee_to_foot*sin(theta_k_swing) geometry.knee_to_foot*cos(theta_k_swing)+geometry.hip_to_knee]) * [hip_x - swing_foot_coords[1]; hip_y - swing_foot_coords[2]]
                
                theta_h_swing = atan(v[1], v[2]) - geometry.theta_hip

                # theta_h_stance = atan(hip_x, hip_y) - atan(-geometry.hip_to_knee-geometry.knee_to_foot*cos(theta_k_stance), geometry.knee_to_foot*sin(theta_k_stance)) - geometry.theta_hip

                # theta_h_swing = atan(hip_x - swing_foot_coords[1], hip_y - swing_foot_coords[2]) - atan(-geometry.hip_to_knee-geometry.knee_to_foot*cos(theta_k_swing), geometry.knee_to_foot*sin(theta_k_swing)) - geometry.theta_hip
                
                config = grounded_left_foot ? SVector(theta_h_stance, theta_k_stance, theta_h_swing, theta_k_swing) : SVector(theta_h_swing, theta_k_swing, theta_h_stance, theta_k_stance)
                rec = MP.get_elem_by_coord(state_grid, config)
                inflated_rec = UT.HyperRectangle(
                    rec.lb - 0.1*(rec.ub - rec.lb),
                    rec.ub + 0.1*(rec.ub - rec.lb)
                )
                push!(target_sets, rec)
            end
        end
    end

    return UT.LazySetUnion(target_sets)

end



end