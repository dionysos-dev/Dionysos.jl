module PathPlanningSgMPC

using StaticArrays
using MathematicalSystems, HybridSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

@enum ApproxMode GROWTH LINEARIZED

function dynamicofsystem()
	# System eq x' = F_sys(x, u)
	function F_sys(x, u)
		return SVector{3}(
			x[1] + u[1] * cos(x[3]),
			x[2] + u[1] * sin(x[3]),
			x[3] + u[2] % (2 * π),
		)
	end

	# We define the growth bound function of $f$:
	ngrowthbound = 5

	# We define the growth bound function of $f$:
	# FIXME: Double check if this is need and is the correct growth bound
	function L_growthbound(u)
		beta = abs(u[1]) # We assume that the velocity is bounded by 1?
		gamma = abs(u[2])  # We assume that the angle is bounded by 2π
		return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, beta, beta, gamma)
	end

	# We define the linearized system of $f$:
	# DF(x, u) = ∂f/∂x(x, u)
	function DF_sys(x, u)
		return SMatrix{3, 3}(
			1.0,
			0.0,
			-u[1] * cos(x[3]),
			0.0,
			1.0,
			-u[1] * sin(x[3]),
			0.0,
			0.0,
			1,
		)
	end

	# We define the bound of the linearized system of $f$:
	bound_DF(u) = 0.0
	bound_DDF(u) = 0.0

	return F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF
end

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


function get_obstacles(
	_X_;
	lb_x1 = -3.5,
	ub_x1 = 3.5,
	lb_x2 = -2.6,
	ub_x2 = 2.6,
	h = 0.1,
)
	# Define the obstacles
	x1 = range(lb_x1, stop = ub_x1, step = h)
	x2 = range(lb_x2, stop = ub_x2, step = h)
	steps1, steps2 = length(x1), length(x2)

	X1 = x1' .* ones(steps2)
	X2 = ones(steps1)' .* x2

	Z1 = (X1 .^ 2 .- X2 .^ 2) .<= 4
	Z2 = (4 .* X2 .^ 2 .- X1 .^ 2) .<= 16


	# Find the upper and lower bounds of X1 and X2 for the obstacle 
	grid = Z1 .& Z2

	return [
		UT.HyperRectangle(SVector(x1[x1lb], x2[x2lb], _X_.lb[3]), SVector(x1[x1ub], x2[x2ub], _X_.ub[3]))
		for (x1lb, x2lb, x1ub, x2ub) in extract_rectangles(grid)
	]

end

function system(
	_X_;
	_U_ = UT.HyperRectangle(SVector(0.2, -1.0), SVector(2.0, 1.0)),
	xdim = 3,
	udim = 2,
	sysnoise = SVector(0.0, 0.0, 0.0),
	measnoise = SVector(0.0, 0.0, 0.0),
	tstep = 0.3,
	nsys = 5,
	approx_mode::ApproxMode = GROWTH,
)
	F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF = dynamicofsystem()
	contsys = nothing
	if approx_mode == GROWTH
		contsys = ST.NewControlSystemGrowthRK4(
			tstep,
			F_sys,
			L_growthbound,
			sysnoise,
			measnoise,
			nsys,
			ngrowthbound,
		)
	elseif approx_mode == LINEARIZED
		contsys = ST.NewControlSystemLinearizedRK4(
			tstep,
			F_sys,
			DF_sys,
			bound_DF,
			bound_DDF,
			measnoise,
			nsys,
		)
	end
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
with the system, the initial and target domains, and null cost functions.
"""
function problem(; simple = false, approx_mode::ApproxMode = GROWTH)
	if simple
		_X_ = UT.HyperRectangle(SVector(-3.5, -2.6, -pi), SVector(3.5, 2.6, pi))
		_I_ = UT.HyperRectangle(SVector(1.0, -1.5, 0.0), SVector(1.2, -1.7, 0.0))
		_T_ = UT.HyperRectangle(SVector(0.5, 0.5, 0.0), SVector(0.7, 0.7, 0.0))
	else
		_X_ = UT.HyperRectangle(SVector(-3.5, -2.6, -pi), SVector(3.5, 2.6, pi))
		_I_ = UT.HyperRectangle(SVector(1.0, -1.5, 0.0), SVector(1.2, -1.7, 0.0))
		_T_ = UT.HyperRectangle(SVector(sqrt(32.0 / 3.0), sqrt(20.0 / 3.0), 0.0), SVector(sqrt(32.0 / 3.0)-0.2, sqrt(20.0 / 3.0)-0.2, 0.0))
	end
	obs = get_obstacles(_X_)
	obstacles_LU = filter_obstacles(_X_, _I_, _T_, obs)
	_X_ = UT.LazySetMinus(_X_, obstacles_LU)
	sys = system(_X_; approx_mode = approx_mode)
	problem = PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
	return problem
end

end
