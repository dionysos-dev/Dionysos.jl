# Tests for robot_problem.jl: the symbolic biped, its kinematics, its two 5-D design models and the
# footstep scenario. Run with an environment that has Dionysos, Symbolics, IntervalArithmetic,
# Plots, MathematicalSystems, LazySets, StaticArrays and MathOptInterface:
#   julia --project=<env> tests/test_robot_problem.jl
using Test
using Random
using LinearAlgebra
using StaticArrays
using Symbolics
using IntervalArithmetic          # loads Dionysos' interval/Symbolics extension
import Dionysos
import LazySets
import MathOptInterface as MOI

const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MP = Dionysos.Mapping

include(joinpath(@__DIR__, "..", "robot_problem.jl"))
include(joinpath(@__DIR__, "..", "LocalGrowthAbstraction.jl"))
import .UnderactuatedRobotProblem as RP
import .LocalGrowthAbstraction as LGA
const MB = RP.Multibody2D

Random.seed!(1234)

# Reference values of the ankle acceleration `q̈₁` of the velocity-controlled design model
# `x = [q̇₁, q₁, q₂, q₃, q₄]`, `u = [v₂, v₃, v₄]` (counter-clockwise angles), computed once with
# RigidBodyDynamics on `four_link_biped.urdf` (angles converted to counter-clockwise). The URDF and
# the symbolic model describe the same robot, so the dynamics must reproduce them.
const RBD_REFERENCE = [
    ([0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0),
    ([0.5, -0.1, 0.2, -0.5, 0.0], [1.0, -1.0, 1.0], 2.11277496402945),
    ([-1.2, 0.3, -0.4, 0.8, -0.6], [2.0, 0.0, -2.0], 0.379467317406006),
    ([0.0, -0.05, 0.0, -0.2, -0.5], [0.0, 2.0, 0.0], -0.125279250230014),
    ([1.5, 0.2, 0.7, 0.4, -1.1], [-1.5, 1.2, 0.3], 11.7454250700264),
    ([-0.3, -0.4, -0.8, 1.0, -0.2], [0.7, -2.0, 2.0], -25.2365818729479),
]

# Closed-form forward kinematics (counter-clockwise absolute angles from the upward vertical),
# derived by hand: independent of the symbolic model.
function closed_form(q, L1, L2)
    φ1 = q[1]
    φ2 = φ1 + q[2]
    φ3 = φ2 + q[3]
    φ4 = φ3 + q[4]
    knee = L2 * [-sin(φ1), cos(φ1)]
    hip = knee + L1 * [-sin(φ2), cos(φ2)]
    swing_knee = hip - L1 * [-sin(φ3), cos(φ3)]
    swing_foot = swing_knee - L2 * [-sin(φ4), cos(φ4)]
    return [[0.0, 0.0], knee, hip, swing_knee, swing_foot]
end

const MODEL = RP.symbolic_model()
const VALUES = RP.default_parameters()

# Numeric first row of the mass matrix, `q -> D₁ⱼ(q)`.
const D_ROW1 = first(
    Symbolics.build_function(
        RP.substitute_parameters(MODEL, MB.mass_matrix(MODEL.sys)[1, :], VALUES),
        MODEL.sys.q;
        expression = Val(false),
    ),
)

@testset "UnderactuatedRobotProblem" begin
    @testset "velocity-controlled dynamics reproduce RigidBodyDynamics on the URDF" begin
        dynamics = RP.VelocityControlledDynamics(; model = MODEL)
        for (x, u, expected) in RBD_REFERENCE
            ẋ = dynamics(SVector{5}(x), SVector{3}(u))
            @test ẋ[1] ≈ expected rtol = 1e-9 atol = 1e-12
            @test ẋ[2:5] == SVector(x[1], u...)                 # (q̇₁, v₂, v₃, v₄)
            @test dynamics(x, u)[1] ≈ expected rtol = 1e-9 atol = 1e-12     # plain vectors too
        end
    end

    @testset "model structure and parameters" begin
        @test keys(MODEL.parameters) == keys(VALUES) == (:L1, :L2, :m1, :m2, :mM, :g)
        @test MODEL.nf.actuated == [2, 3, 4] && MODEL.nf.unactuated == [1]
        @test length(MODEL.nf.state) == 8 && length(MODEL.nf.input) == 3
        # The numbers are an argument: another motor mass gives another acceleration.
        x, u = [0.5, -0.1, 0.2, -0.5, 0.0], [1.0, -1.0, 1.0]
        heavy = merge(VALUES, (; mM = 4.0))
        @test RP.VelocityControlledDynamics(; model = MODEL, values = heavy)(x, u)[1] !=
              RP.VelocityControlledDynamics(; model = MODEL)(x, u)[1]
        # `substitute_parameters` leaves nothing but the coordinates and the inputs.
        f = RP.substitute_parameters(MODEL, MB.mass_matrix(MODEL.sys)[1, 1], VALUES)
        @test all(v -> any(w -> isequal(v, w), MODEL.sys.q), Symbolics.get_variables(f))
    end

    @testset "kinematics agree with the closed form, for any lengths" begin
        for values in (VALUES, merge(VALUES, (; L1 = 0.31, L2 = 0.12)))
            geometry = RP.Geometry(MODEL, values)
            @test (geometry.femur, geometry.tibia) == (values.L1, values.L2)
            for k in 1:5
                q = SVector{4}(0.7 * sin(12.9898 * k + 78.233 * j) for j in 1:4)
                points = RP.cartesian_coordinates(q, geometry)
                ref = closed_form(q, values.L1, values.L2)
                @test length(points) == 5
                for (p, r) in zip(points, ref)
                    @test collect(p) ≈ r atol = 1e-12
                end
                @test RP.swing_foot_position(q, geometry) ≈ ref[5] atol = 1e-12
            end
        end
    end

    @testset "static equilibrium has no ankle acceleration" begin
        # With the ankle at rest the gravity torque vanishes when the centre of mass is above the
        # ankle: bisect q₁ for a posture with the swing leg behind, then check q̈₁ = 0 there.
        dynamics = RP.VelocityControlledDynamics(; model = MODEL)
        f(x) = dynamics(x, zeros(3))[1]
        posture(q1) = [0.0, q1, 0.0, -0.2 - q1, 0.0]
        a, b = -0.5, 0.5
        sign(f(posture(a))) != sign(f(posture(b))) || error("no equilibrium")
        for _ in 1:60
            m = (a + b) / 2
            sign(f(posture(m))) == sign(f(posture(a))) ? (a = m) : (b = m)
        end
        @test abs(f(posture((a + b) / 2))) < 1e-9
    end

    @testset "normal-form model and velocity-controlled model are the same system" begin
        reduced = RP.reduced_dynamics(MODEL)
        # State (p₁, q), input v_a, and no dependence on the torques.
        @test isequal(reduced.x[1], MODEL.nf.state[8]) &&
              isequal(reduced.x[2:5], MODEL.nf.state[1:4])
        @test isequal(reduced.u, MODEL.nf.state[5:7])
        @test all(
            e -> isempty(intersect(Symbolics.get_variables(e), MODEL.nf.input)),
            reduced.f,
        )

        f_nf = RP.substitute_parameters(MODEL, reduced.f, VALUES)
        nf = LGA.SymbolicDynamics(f_nf, reduced.x, reduced.u)
        vel = RP.VelocityControlledDynamics(; model = MODEL)
        δ = 1e-4
        for k in 1:6
            q = [0.4 * sin(1.7k + j) for j in 1:4]
            u = SVector{3}(1.5 * cos(2.3k + j) for j in 1:3)
            x_vel = SVector{5}(0.8 * sin(0.9k), q...)             # (q̇₁, q)
            p1(x) = dot(D_ROW1(x[2:5]), [x[1]; u])               # p₁ = (D q̇)₁
            x_nf = SVector{5}(p1(x_vel), q...)                   # the same state, with p₁
            f = nf.f(x_nf, u)
            # q̇₁ = (p₁ − D₁ₐ v_a) / D₁₁ is the velocity we started from, and the rest is v_a.
            @test f[2:5] ≈ [x_vel[1]; u] rtol = 1e-9 atol = 1e-12
            # ṗ₁ = d/dt (D q̇)₁ along the velocity-controlled flow, by a central difference.
            flow(t) = ST.runge_kutta4(vel, x_vel, u, t, 20)
            ṗ1 = (p1(flow(δ)) - p1(flow(-δ))) / (2δ)
            @test f[1] ≈ ṗ1 rtol = 1e-5 atol = 1e-7
        end
    end

    @testset "local growth abstraction of the normal-form model (tiny grid)" begin
        reduced = RP.reduced_dynamics(MODEL)
        f_nf = RP.substitute_parameters(MODEL, reduced.f, VALUES)
        disc = RP.default_discretization(; dx = 0.1, dq1 = 0.1, dq̇1 = 0.05, tstep = 0.05)
        h = MP.get_h(disc.state_grid)
        X = LazySets.Hyperrectangle(; low = -1.5 .* h, high = 1.5 .* h)   # 3 cells per axis
        U = LazySets.Hyperrectangle(;
            low = fill(-disc.u_max, 3),
            high = fill(disc.u_max, 3),
        )
        function build(; kwargs...)
            lgb = LGA.LocalGrowthBound(
                f_nf,
                reduced.x,
                reduced.u,
                X,
                U,
                disc.tstep;
                kwargs...,
            )
            optimizer = LGA.build_abstraction(
                lgb;
                state_grid = disc.state_grid,
                input_grid = disc.input_grid,
                print_level = 0,
                execution_backend = SY.SequentialBackend(),
            )
            return lgb, MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
        end
        lgb, sym = build()
        _, sym_global = build(; refinements = 0)       # the global growth bound, for comparison
        @test SY.get_n_state(sym) == 3^5 && SY.get_n_input(sym) == 27
        function successor_sets(abstract_system)
            sets = Dict{Tuple{Int, Int}, Set{Int}}()
            for (t, s, i) in SY.enum_transitions(abstract_system)
                push!(get!(sets, (s, i), Set{Int}()), t)
            end
            return sets
        end
        successors, successors_global = successor_sets(sym), successor_sets(sym_global)
        @test !isempty(successors) && !isempty(successors_global)
        # The local image of a (cell, input) pair is inside its global image, so the local
        # abstraction enables every pair the global one enables, with fewer successors. (It can
        # have more transitions in total: pairs whose global image leaves the domain are enabled.)
        @test all(
            k -> haskey(successors, k) && successors[k] ⊆ successors_global[k],
            keys(successors_global),
        )

        # Soundness by sampling: a random point of a random enabled cell, simulated over one
        # sampling period, lands in one of the successors of that (cell, input) pair.
        pairs_ = collect(keys(successors))
        violations = 0
        for _ in 1:300
            (s, i) = rand(pairs_)
            x0 = SY.get_concrete_state(sym, s)
            u = SY.get_concrete_input(sym, i)
            p = x0 .+ (2 .* rand(5) .- 1) .* h ./ 2
            y = ST.runge_kutta4(lgb.dynamics.f, SVector{5}(p), u, lgb.tstep, 200)
            t = SY.get_abstract_state(sym, y)
            (t === nothing || t in successors[(s, i)]) || (violations += 1)
        end
        @test violations == 0
    end

    @testset "footstep scenario: carving, target, step problem" begin
        geometry = RP.Geometry(MODEL, VALUES)
        triangle =
            LazySets.VPolygon([SVector(-0.04, 0.0), SVector(0.0, 0.03), SVector(0.04, 0.0)])
        disc = RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1)
        lo = SVector(-0.3, -0.3, -0.5, -0.8, -1.0)
        hi = SVector(0.3, 0.2, 0.5, 1.0, 0.2)
        X_box = LazySets.Hyperrectangle(; low = lo, high = hi)
        removed = RP.infeasible_cells(X_box, disc.state_grid, [triangle], geometry)
        @test !isempty(removed)

        # The carving is sound: no configuration of a kept cell puts the swing foot in the obstacle.
        # The abstraction keeps the cells lying entirely inside the box (`MP.INNER`), as does the
        # carving, so cells straddling its boundary are not states and are not sampled.
        h = MP.get_h(disc.state_grid)
        inside(pos) = (
            c = MP.get_coord_by_pos(disc.state_grid, pos);
            all(lo .- 1e-9 .<= c .- h ./ 2) && all(c .+ h ./ 2 .<= hi .+ 1e-9)
        )
        kept_samples = kept_hits = 0
        for _ in 1:4000
            x = lo .+ rand(5) .* (hi .- lo)
            pos = MP.get_pos_by_coord(disc.state_grid, x)
            (inside(pos) && !(pos in removed)) || continue
            kept_samples += 1
            kept_hits += RP.swing_foot_position(SVector{4}(x[2:5]), geometry) ∈ triangle
        end
        @test kept_samples > 100 && kept_hits == 0
        # The ankle coordinate does not matter: the removed cells are whole columns along it.
        pos = first(removed.positions)
        @test all((k, pos[2:5]...) in removed.positions for k in -3:3)

        # The inverse kinematics place the swing foot on the foothold, with the knee bent backwards.
        foothold = SVector(0.2, 0.0)
        reachable = 0
        for q1 in -0.3:0.1:0.2, q2 in -0.5:0.1:0.5
            hip = RP.cartesian_coordinates(SVector(q1, q2, 0.0, 0.0), geometry)[3]
            ik = RP._swing_leg_ik(foothold - hip, geometry)
            ik === nothing && continue                # the foothold is out of reach of the swing leg
            φ3, φ4 = ik
            q = SVector(q1, q2, φ3 - (q1 + q2), φ4 - φ3)
            @test RP.swing_foot_position(q, geometry) ≈ foothold atol = 1e-12
            @test q[4] <= 0
            reachable += 1
        end
        @test reachable > 5

        # Target cells: the ankle coordinate is free by default, or restricted to whole cells.
        all_cells = RP.target_cells(disc.state_grid, foothold, geometry, X_box)
        ankle_all = Set(c[1] for c in all_cells)
        @test !isempty(all_cells) && issubset(Set(-2:2), ankle_all)     # every inner ankle cell
        some = RP.target_cells(
            disc.state_grid,
            foothold,
            geometry,
            X_box;
            velocity = (-0.15, 0.15),
        )
        @test Set(c[1] for c in some) == Set(-1:1)
        @test_throws ErrorException RP.target_cells(
            disc.state_grid,
            foothold,
            geometry,
            X_box;
            velocity = (0.01, 0.02),
        )

        domain = RP.RobotDomainConfig(;
            x_lb = lo,
            x_ub = hi,
            u_lb = SVector(-1.0, -1.0, -1.0),
            u_ub = SVector(1.0, 1.0, 1.0),
        )
        system = RP.system(;
            dynamics = RP.VelocityControlledDynamics(; model = MODEL),
            domain,
            removed_cells = removed,
        )
        q0 = SVector(-0.05, 0.0, -0.2, -0.5)         # swing foot 15 cm behind the stance foot
        problem = RP.step_problem(system, disc.state_grid; q0, foothold, geometry)
        v_init = LazySets.low(problem.initial_set)[1], LazySets.high(problem.initial_set)[1]
        @test v_init == (lo[1], hi[1])                # whole interval by default
        problem = RP.step_problem(
            system,
            disc.state_grid;
            q0,
            foothold,
            geometry,
            initial_velocity = (-0.1, 0.1),
        )
        @test (
            LazySets.low(problem.initial_set)[1],
            LazySets.high(problem.initial_set)[1],
        ) == (-0.1, 0.1)
        # The cell of a start configuration with the swing foot inside the obstacle is carved out.
        @test_throws ErrorException RP.step_problem(
            system,
            disc.state_grid;
            q0 = SVector(0.0, 0.0, 0.0, 0.0),
            foothold,
            geometry,
        )
        @test_throws ErrorException RP.step_problem(
            system,
            disc.state_grid;
            q0,
            foothold,
            geometry,
            initial_velocity = (5.0, 6.0),
        )

        # The figures of a closed-loop run, on a short synthetic trajectory.
        xs = [SVector(0.0, -0.05 + 0.01k, 0.0, -0.2 + 0.05k, -0.5 + 0.1k) for k in 0:6]
        us = [SVector(0.0, 0.5, -0.5) for _ in 1:6]
        mktempdir() do dir
            files = RP.plot_footstep(
                ST.Trajectory(xs; inputs = us),
                disc.state_grid,
                removed;
                x_lb = lo,
                x_ub = hi,
                tstep = disc.tstep,
                geometry,
                obstacle = [triangle],
                foothold,
                file_prefix = joinpath(dir, "step"),
                ankle_label = "q̇₁ (rad/s)",
            )
            @test keys(files) == (:snapshot, :ankle, :animation)
            @test all(f -> isfile(f) && filesize(f) > 0, values(files))
            # A trajectory without action has nothing to draw.
            @test RP.plot_footstep(
                ST.Trajectory(xs[1:1]; inputs = SVector{3, Float64}[]),
                disc.state_grid,
                removed;
                x_lb = lo,
                x_ub = hi,
                tstep = disc.tstep,
                file_prefix = joinpath(dir, "none"),
            ) === nothing
            @test !isfile(joinpath(dir, "none.png"))
        end
    end
end
