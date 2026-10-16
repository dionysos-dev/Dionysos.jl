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
import MathematicalSystems
import MathOptInterface as MOI

const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MP = Dionysos.Mapping
const PR = Dionysos.Problem
const UT = Dionysos.Utils

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

# The velocity-controlled model written with the ankle *velocity*, x = [q̇₁, q₁…q₄], u = [v₂, v₃, v₄]:
# ẋ = [q̈₁, q̇₁, v₂, v₃, v₄] with q̈₁ = h₁ / D₁₁, the first row of D q̈ = h + Bτ (whose row of B is
# zero) with zero actuated accelerations. It is derived from the mass matrix and the bias vector
# alone, independently of the normal form of `RP.reduced_dynamics`, and is what the RigidBodyDynamics
# reference values are expressed in.
function velocity_model(model = MODEL, values = VALUES)
    sys = model.sys
    q̈1 = RP.substitute_parameters(
        model,
        MB.bias_vector(sys)[1] / MB.mass_matrix(sys)[1, 1],
        values,
    )
    f = Symbolics.build_function(
        q̈1,
        [sys.qd[1]; sys.q],
        sys.qd[2:4];
        expression = Val(false),
        cse = true,
    )
    q̈1_function = f isa Tuple ? first(f) : f
    return (x, u) -> SVector(q̈1_function(x, u), x[1], u[1], u[2], u[3])
end

@testset "UnderactuatedRobotProblem" begin
    @testset "velocity-controlled dynamics reproduce RigidBodyDynamics on the URDF" begin
        dynamics = velocity_model()
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
        @test velocity_model(MODEL, heavy)(x, u)[1] != velocity_model()(x, u)[1]
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
        dynamics = velocity_model()
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
        vel = velocity_model()
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

    @testset "input set: integer multiples of du up to input_steps · du" begin
        d1 = RP.default_discretization(; dx = 0.1, dq1 = 0.1, dq̇1 = 0.05, tstep = 0.05)
        @test d1.input_steps == 1 && d1.u_max == d1.du == 2.0
        @test RP.input_count(d1) == 27 && RP.input_count(d1; one_joint_per_step = true) == 7
        @test_throws ErrorException RP.default_discretization(; input_steps = 0)
        for k in (2, 3)
            dk = RP.default_discretization(;
                dx = 0.1,
                dq1 = 0.1,
                dq̇1 = 0.05,
                tstep = 0.05,
                input_steps = k,
            )
            @test dk.u_max == k * dk.du && MP.get_h(dk.input_grid) == fill(dk.du, 3)
            @test RP.input_count(dk) == (2k + 1)^3
            @test RP.input_count(dk; one_joint_per_step = true) == 1 + 6k
        end
        # The abstraction really has these inputs, all multiples of du: 5³ for input_steps = 2.
        reduced = RP.reduced_dynamics(MODEL)
        f_nf = RP.substitute_parameters(MODEL, reduced.f, VALUES)
        disc = RP.default_discretization(;
            dx = 0.1,
            dq1 = 0.1,
            dq̇1 = 0.05,
            tstep = 0.05,
            input_steps = 2,
        )
        h = MP.get_h(disc.state_grid)
        X = LazySets.Hyperrectangle(; low = -0.5 .* h, high = 0.5 .* h)   # one cell
        U = LazySets.Hyperrectangle(;
            low = fill(-disc.u_max, 3),
            high = fill(disc.u_max, 3),
        )
        lgb = LGA.LocalGrowthBound(f_nf, reduced.x, reduced.u, X, U, disc.tstep)
        optimizer = LGA.build_abstraction(
            lgb;
            state_grid = disc.state_grid,
            input_grid = disc.input_grid,
            execution_backend = SY.SequentialBackend(),
        )
        sym = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
        @test SY.get_n_input(sym) == RP.input_count(disc) == 125
        inputs = [SY.get_concrete_input(sym, i) for i in 1:SY.get_n_input(sym)]
        @test all(
            u -> all(c -> isapprox(c / disc.du, round(c / disc.du); atol = 1e-9), u),
            inputs,
        )
        @test extrema(maximum(abs, u) for u in inputs) == (0.0, disc.u_max)
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
            return lgb,
            optimizer,
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
        end
        lgb, optimizer, sym = build()
        _, _, sym_global = build(; refinements = 0)       # the global growth bound, for comparison
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

        # Start state and closed loop, on the same abstraction: from the origin, reach the cells one
        # step up along the grounded knee `q₂` with any other coordinate.
        origin_q = SVector(0.0, 0.0, 0.0, 0.0)
        x0 = RP.start_state(disc.state_grid, origin_q, X)
        @test x0 == zeros(5)
        # The first coordinate is the middle of the range; the configuration is snapped to its cell.
        x1 = RP.start_state(
            disc.state_grid,
            SVector(0.04, 0.0, 0.0, 0.0),
            X;
            velocity = (0.04, 0.06),
        )
        @test x1 ≈ [0.05, 0.0, 0.0, 0.0, 0.0]
        # A start in a carved-out cell is refused; a carved cell elsewhere is not an issue.
        @test_throws ErrorException RP.start_state(
            disc.state_grid,
            origin_q,
            X;
            removed = MP.CellUnion(disc.state_grid, [(0, 0, 0, 0, 0)]),
        )
        @test x0 == RP.start_state(
            disc.state_grid,
            origin_q,
            X;
            removed = MP.CellUnion(disc.state_grid, [(1, 1, 1, 1, 1)]),
        )

        target = MP.CellUnion(
            disc.state_grid,
            [(a, i, 1, j, k) for a in -1:1, i in -1:1, j in -1:1, k in -1:1],
        )
        initial = LazySets.Hyperrectangle(; low = x0 .- 0.25 .* h, high = x0 .+ 0.25 .* h)
        problem = PR.OptimalControlProblem(
            ST.get_system(lgb),
            initial,
            target,
            nothing,
            nothing,
            PR.Infinity(),
        )
        MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
        MOI.optimize!(optimizer)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
        traj = redirect_stdout(devnull) do
            return RP.closed_loop(optimizer, problem, x0; horizon = 10, ankle_name = "p₁")
        end
        @test length(ST.inputs(traj)) == 1 && last(ST.states(traj)) ∈ target
    end

    @testset "cell counts of the boxes of the experiment" begin
        grid =
            RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05).state_grid
        box(lo, hi) = LazySets.Hyperrectangle(; low = SVector(lo...), high = SVector(hi...))
        # The default box of the experiment: 19 × 9 × 9 × 17 × 11 cells.
        narrow = box((-0.2, -0.3, -0.5, -0.8, -1.0), (0.2, 0.2, 0.5, 1.0, 0.2))
        @test RP.abstraction_size(narrow, grid).cells_in_box ==
              19 * 9 * 9 * 17 * 11 ==
              287_793
        # A box with ±0.5 on the ankle, ±0.5 on q₁ and ±1 on the swing joints looks only "slightly"
        # wider but has ten times more cells (and 17 times more pairs).
        wide = box((-0.5, -0.5, -0.5, -1.0, -1.0), (0.5, 0.5, 0.5, 1.0, 1.0))
        @test RP.abstraction_size(wide, grid).cells_in_box ==
              49 * 19 * 9 * 19 * 19 ==
              3_024_819
        # Cells straddling the boundary are not kept: a box not aligned with the grid loses them.
        @test RP.abstraction_size(box(fill(-0.06, 5), fill(0.06, 5)), grid).cells_in_box ==
              5 * 1 * 1 * 1 * 1
    end

    @testset "footstep scenario: carving, target, step problem" begin
        geometry = RP.Geometry(MODEL, VALUES)
        triangle =
            LazySets.VPolygon([SVector(-0.04, 0.0), SVector(0.0, 0.03), SVector(0.04, 0.0)])
        disc = RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1)
        lo = SVector(-0.3, -0.3, -0.5, -0.8, -1.0)
        hi = SVector(0.3, 0.2, 0.5, 1.0, 0.2)
        X_box = LazySets.Hyperrectangle(; low = lo, high = hi)
        carving = RP.infeasible_cells(
            X_box,
            disc.state_grid,
            [triangle],
            geometry;
            obstacle_name = "lift constraint",
        )
        removed = carving.cells
        @test !isempty(removed)

        # The size of the abstraction, before computing it: the cells lying entirely inside the box
        # (what the abstraction keeps) minus the carved ones, as counted by brute force.
        inner = collect(MP.get_pos_from_set(disc.state_grid, X_box, MP.INNER))
        sized = RP.abstraction_size(X_box, disc.state_grid, carving; threads = 10)
        @test sized.cells_in_box == length(inner)
        @test sized.cells == count(pos -> !(pos in removed), inner)
        @test RP.abstraction_size(X_box, disc.state_grid, removed).cells == sized.cells
        @test sized.pairs == 27 * sized.cells
        @test sized.transitions == RP.TRANSITIONS_PER_PAIR * sized.pairs
        @test sized.minutes ≈ sized.pairs * RP.CPU_SECONDS_PER_PAIR / 10 / 60
        @test RP.abstraction_size(X_box, disc.state_grid; ninputs = 7).pairs ==
              7 * length(inner)

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

        # Which constraint removed what: the ground and the named obstacle, and every removed cell
        # is accounted for by at least one of them.
        @test first.(carving.causes) == ["ground", "lift constraint"]
        ground, lift = last.(carving.causes)
        @test !isempty(lift)
        @test isempty(ground)       # no posture of this box puts the foot surely below the ground
        @test all(
            !isempty(RP.carved_by(carving, p)) for
            p in Iterators.take(removed.positions, 500)
        )
        @test RP.carved_by(carving, (0, first(lift)...)) ⊇ ["lift constraint"]
        @test RP.carved_by(carving, (1, 99, 99, 99, 99)) == String[]
        # The default name, one name per obstacle, and the error on a wrong number of names.
        far = LazySets.VPolygon([SVector(0.3, 0.0), SVector(0.35, 0.03), SVector(0.4, 0.0)])
        @test first.(
            RP.infeasible_cells(X_box, disc.state_grid, triangle, geometry).causes,
        ) == ["ground", "obstacle"]
        @test first.(
            RP.infeasible_cells(X_box, disc.state_grid, nothing, geometry).causes,
        ) == ["ground"]
        two = RP.infeasible_cells(
            X_box,
            disc.state_grid,
            [triangle, far],
            geometry;
            obstacle_name = ["lift constraint", "step"],
        )
        @test first.(two.causes) == ["ground", "lift constraint", "step"]
        @test last.(two.causes)[2] == lift                # the same obstacle, the same cells
        @test two.cells.positions ⊇ removed.positions
        @test_throws ErrorException RP.infeasible_cells(
            X_box,
            disc.state_grid,
            [triangle, far],
            geometry;
            obstacle_name = ["lift constraint"],
        )

        # A start in a carved-out cell is refused, and the error names the constraint that did it.
        message(f) =
            try
                f()
                ""
            catch e
                e isa ErrorException ? e.msg : rethrow()
            end
        start_in(pos4; removed = carving) = message() do
            return RP.start_state(
                disc.state_grid,
                SVector(MP.get_coord_by_pos(disc.state_grid, (0, pos4...))[2:5]...),
                X_box;
                removed,
                geometry,
            )
        end
        only_lift = first(lift)
        msg = start_in(only_lift)
        @test occursin("carved out by \"lift constraint\"", msg) && !occursin("ground", msg)
        @test occursin("within 7.5 cm", msg) && occursin("(x, z) = (", msg)
        # The ground constraint, in a box of strongly bent postures (the stance knee bent by 1 rad with
        # the swing femur vertical: the foot is 9 cm below the ground), on a finer grid.
        fine = RP.default_discretization(; dx = 0.05, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.05)
        X_low = LazySets.Hyperrectangle(;
            low = SVector(-0.1, -0.1, 0.8, -1.0, -0.1),
            high = SVector(0.1, 0.1, 1.0, -0.8, 0.1),
        )
        low = RP.infeasible_cells(X_low, fine.state_grid, [triangle], geometry)
        only_ground = first(last(low.causes[1]))
        @test low.causes[1].first == "ground" && !isempty(last(low.causes[1]))
        msg = message() do
            return RP.start_state(
                fine.state_grid,
                SVector(MP.get_coord_by_pos(fine.state_grid, (0, only_ground...))[2:5]...),
                X_low;
                removed = low,
                geometry,
            )
        end
        @test occursin("carved out by \"ground\"", msg) && !occursin("obstacle", msg)
        @test occursin("below the ground", msg)
        # A cell removed by two constraints names both.
        both = RP.Carving(
            removed,
            ["ground" => Set([only_lift]), "lift constraint" => Set([only_lift])],
        )
        @test occursin("\"ground\" (", start_in(only_lift; removed = both)) &&
              occursin(" and \"lift constraint\" (", start_in(only_lift; removed = both))
        # Only the cells are known: the generic message, still an error.
        @test occursin("carved out", start_in(only_lift; removed = removed))
        @test !occursin("lift constraint", start_in(only_lift; removed = removed))
        # A cell that is not carved out is fine.
        @test RP.start_state(
            disc.state_grid,
            SVector(-0.05, 0.0, -0.2, -0.5),
            X_box;
            removed = carving,
            geometry,
        ) isa SVector{5}

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
        # A foothold tolerance makes the target thicker: the cells for the exact foothold, and more.
        thick =
            RP.target_cells(disc.state_grid, foothold, geometry, X_box; tolerance = 0.03)
        @test Set(all_cells) ⊆ Set(thick) && length(thick) > 2 * length(all_cells)
        thicker =
            RP.target_cells(disc.state_grid, foothold, geometry, X_box; tolerance = 0.06)
        @test Set(thick) ⊆ Set(thicker) && length(thicker) > length(thick)
        @test_throws ErrorException RP.target_cells(
            disc.state_grid,
            foothold,
            geometry,
            X_box;
            tolerance = -0.01,
        )
        @test_throws ErrorException RP.target_cells(
            disc.state_grid,
            foothold,
            geometry,
            X_box;
            velocity = (0.01, 0.02),
        )

        # `step_problem` only reads the domain of the system: a placeholder with the carved domain.
        system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
            (x, u) -> x,
            5,
            3,
            UT.set_minus(X_box, removed),
            LazySets.Hyperrectangle(; low = fill(-1.0, 3), high = fill(1.0, 3)),
        )
        q0 = SVector(-0.05, 0.0, -0.2, -0.5)         # swing foot 15 cm behind the stance foot
        problem = RP.step_problem(system, disc.state_grid; q0, foothold, geometry)
        v_init = LazySets.low(problem.initial_set)[1], LazySets.high(problem.initial_set)[1]
        @test v_init == (lo[1], hi[1])                # whole interval by default
        tolerant = RP.step_problem(
            system,
            disc.state_grid;
            q0,
            foothold,
            geometry,
            foothold_tolerance = 0.03,
        )
        @test length(tolerant.target_set.positions) > length(problem.target_set.positions)
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
            # With the dynamics, every step is re-integrated at intermediate instants (the same
            # figures, from a refined trajectory). Constant motor velocities, a still ankle.
            refined = RP.plot_footstep(
                ST.Trajectory(xs; inputs = us),
                disc.state_grid,
                removed;
                x_lb = lo,
                x_ub = hi,
                tstep = disc.tstep,
                geometry,
                obstacle = [triangle],
                foothold,
                file_prefix = joinpath(dir, "refined"),
                ankle_label = "q̇₁ (rad/s)",
                dynamics = (x, u) -> SVector(0.0, 0.0, u[1], u[2], u[3]),
                subframes = 4,
            )
            @test all(f -> isfile(f) && filesize(f) > 0, values(refined))
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
