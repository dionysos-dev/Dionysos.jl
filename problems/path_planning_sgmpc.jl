module PathPlanningSgMPC

using StaticArrays, FillArrays
using MathematicalSystems, HybridSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

# TODO: Move to Dionysos.Utils (the function is already there, but it imposes that Q, R, and N must be the same type and *size* and that q and r should be the same type and *size*)
struct QuadraticStateControlFunction{T} <: UT.ScalarControlFunction
    Q::AbstractMatrix{T}
    R::AbstractMatrix{T}
    N::AbstractMatrix{T}
    q::AbstractArray{T}
    r::AbstractArray{T}
    v::T

    #=
    function QuadraticStateControlFunction(Q::MT, R::MT, N::MT, q::AT, r::AT, v::T) where {T, MT <: AbstractMatrix{T}, AT <: AbstractArray{T}}
       # Perform size checks
       @assert size(Q, 1) == size(Q, 2) "Q must be square"
       @assert size(R, 1) == size(R, 2) "R must be square"
       @assert size(Q, 1) == size(N, 1) "Q and N must have compatible dimensions"
       @assert size(R, 1) == size(N, 2) "R and N must have compatible dimensions"
       @assert size(Q, 1) == length(q) "Q and q must have compatible dimensions"
       @assert size(R, 1) == length(r) "R and r must have compatible dimensions"
       new{T, MT, AT}(Q, R, N, q, r, v)
    end
    =#

end
function function_value(f::QuadraticStateControlFunction, x, u)
    return x'f.Q * x + u'f.R * u + 2 * (x'f.N * u + x'f.q + u'f.r) + f.v
end
function get_full_psd_matrix(f::QuadraticStateControlFunction)
    return [
        f.Q f.N f.q
        f.N' f.R f.r
        f.q' f.r' f.v
    ]
end
##########

function filter_obstacles(_X_, _I_, _T_, obs)
    obstacles = typeof(_X_)[]
    for ob in obs
        if ob ⊆ _X_ && isempty(ob ∩ _I_) && isempty(ob ∩ _T_)
            push!(obstacles, ob)
        end
    end
    obstacles_LU = UT.LazyUnionSetArray(obstacles)
    return obstacles_LU
end

function extract_rectangles(matrix)
    if isempty(matrix)
        return []
    end

    n, m = size(matrix)
    tlx, tly, brx, bry = Int[], Int[], Int[], Int[]

    # Build histogram heights
    for i in 1:n
        j = 1
        while j <= m
            if matrix[i, j] == 1
                j += 1
                continue
            end
            push!(tlx, j)
            push!(tly, i)
            while j <= m && matrix[i, j] == 0
                j += 1
            end
            push!(brx, j - 1)
            push!(bry, i)
        end
    end

    return zip(tlx, tly, brx, bry)
end

function get_obstacles(_X_; lb_x1 = -3.5, ub_x1 = 3.5, lb_x2 = -2.6, ub_x2 = 2.6, h = 0.1)
    # Define the obstacles
    x1 = range(lb_x1; stop = ub_x1, step = h)
    x2 = range(lb_x2; stop = ub_x2, step = h)
    steps1, steps2 = length(x1), length(x2)

    X1 = x1' .* ones(steps2)
    X2 = ones(steps1)' .* x2

    Z1 = (X1 .^ 2 .- X2 .^ 2) .<= 4
    Z2 = (4 .* X2 .^ 2 .- X1 .^ 2) .<= 16

    # Find the upper and lower bounds of X1 and X2 for the obstacle 
    grid = Z1 .& Z2

    return [
        UT.HyperRectangle(
            SVector(x1[x1lb], x2[x2lb], _X_.lb[3]),
            SVector(x1[x1ub], x2[x2ub], _X_.ub[3]),
        ) for (x1lb, x2lb, x1ub, x2ub) in extract_rectangles(grid)
    ]
end

function system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(0.2, -1.0), SVector(2.0, 1.0)),
    xdim = 3,
    udim = 2,
    sysnoise = SVector(0.0, 0.0, 0.0),
    measnoise = SVector(0.0, 0.0, 0.0),
    tstep = 1.0,
    nsys = 5,
    ngrowthbound = 5,
)
    sys_map = let nsys = nsys
        (x, u, _) -> SVector{3}(
            x[1] + u[1] * cos(x[3]),
            x[2] + u[1] * sin(x[3]),
            (x[3] + u[2]) % (2 * π),
        )
    end
    sys_inv_map = let nsys = nsys
        (x, u, _) -> SVector{3}(
            x[1] - u[1] * cos((x[3] - u[2]) % (2 * π)),
            x[2] - u[1] * sin((x[3] - u[2]) % (2 * π)),
            (x[3] - u[2]) % (2 * π), # to check
        )
    end
    growthbound_map = let ngrowthbound = ngrowthbound
        (r, u, _) -> SVector{3}(u[1] * r[3], u[1] * r[3], 0.0)
    end
    contsys = ST.ControlSystemGrowth(
        tstep,
        sysnoise,
        measnoise,
        sys_map,
        growthbound_map,
        sys_inv_map,
    )
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        contsys,
        xdim,
        udim,
        _X_,
        _U_,
    )
end

""""
	problem()

This function create the system with `PathPlanningSgMPC.system`.

Then, we define initial and target domains for the state of the system.

Finally, we instantiate our Reachability Problem as an OptimalControlProblem 
with the system, the initial and target domains, and cost functions.

The optimization problem is set with a prediction horizon `N = 20` and the stage cost:
```math
l(x, u) = 100 |(x_1, x_2)^T - x_r|^2 + |u|^2
```
The terminal cost is:
```math
L(x) = 100 |(x_1, x_2)^T - x_r|^2
```

These costs are defined by the quadratic form `x' * Q * x + u' * R * u + 2 * (x' * N * u + x' * q + u' * r) + v`.
"""
function problem(;
    sgmpc = false,
    initial = SVector(1.0, -1.7, 0.0),
    target = SVector(0.5, 0.5, -pi),
)
    _X_ = UT.HyperRectangle(SVector(-3.5, -2.6, -pi), SVector(3.5, 2.6, pi))
    _I_ = UT.HyperRectangle(initial, initial + SVector(0.2, 0.2, 0.0))
    _T_ = UT.HyperRectangle(target, target + SVector(0.2, 0.2, 2 * pi))

    obs = get_obstacles(_X_)
    obstacles_LU = filter_obstacles(_X_, _I_, _T_, obs)
    _X_ = UT.LazySetMinus(_X_, obstacles_LU)
    sys = system(_X_)

    if sgmpc
        Q = SMatrix{3, 3}(100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0)
        R = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)
        N = zeros(Float64, 3, 2)

        target_without_last_element = SVector{3}(1.0, 1.0, 0.0) .* target
        q = -100 * target_without_last_element
        r = SVector{2}(0.0, 0.0)
        v = 100 * target_without_last_element' * target_without_last_element
        NT = 20

        zero_cost = QuadraticStateControlFunction(
            zeros(Float64, 3, 3),
            zeros(Float64, 2, 2),
            zeros(Float64, 3, 2),
            [0.0, 0.0, 0.0],
            [0.0, 0.0],
            0.0,
        )
        state_cost = [QuadraticStateControlFunction(Q, R, N, q, r, v) for _ in 1:(NT - 1)]
        terminal_cost = [zero_cost for _ in 1:NT]
        terminal_cost[NT] =
            QuadraticStateControlFunction(Q, zeros(Float64, 2, 2), N, q, r, v)

        problem = PB.OptimalControlProblem(sys, _I_, _T_, state_cost, terminal_cost, NT)
    else
        #TODO: Convert LazySetMinus to HyperRectangle to use SafetyProblem
        #problem = PB.SafetyProblem(sys, sys.X, sys.X, PB.Infinity())
        problem = PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    end

    return problem
end

end
