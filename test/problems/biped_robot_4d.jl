module TestBipedRobot4D

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import Random
import MathOptInterface as MOI

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "BipedRobot",
        "4D_model",
        "robot_problem.jl",
    ),
)
import .RobotProblem as RP

const GEOM = RP.default_geometry()
const RNG = Random.MersenneTwister(1)

rand_config(rng; bar = 0.5) = SVector{4}(bar .* (2 .* rand(rng, 4) .- 1))

@testset "FK/IK round trip" begin
    for grounded_left_foot in (true, false)
        for _ in 1:100
            θ = rand_config(RNG)
            cartesian = RP.cartesian_coordinates(θ, GEOM, grounded_left_foot)
            θ_back, g_back = RP.angular_coordinates(cartesian, GEOM)
            @test g_back == grounded_left_foot
            @test θ_back ≈ θ atol = 1e-9
        end
    end
end

@testset "swing-foot Lipschitz bound" begin
    # The deviation bound must dominate the actual foot displacement between
    # any two configurations of a cell.
    r = SVector(0.05, 0.05, 0.05, 0.05)
    dev = RP.swing_foot_deviation_bound(GEOM, r)
    for _ in 1:200
        c = rand_config(RNG)
        δ = SVector{4}(r .* (2 .* rand(RNG, 4) .- 1))
        f_c = RP.swing_foot_position(c, GEOM, true)
        f_x = RP.swing_foot_position(c + δ, GEOM, true)
        @test sqrt(sum((f_x - f_c) .^ 2)) <= dev + 1e-12
    end
end

const DISC = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)
# The obstacle sits *behind* the swing trajectory (which goes from x ≈ -0.15 to
# x = 0.2): sound carving at dx = 0.1 carries a Lipschitz margin of ~5.5 cm, so
# an obstacle placed under the swing path would provably disconnect the free
# space at this resolution — actually climbing over needs a finer grid (see the
# examples/BipedRobot driver). Carving and inter-sample safety are exercised
# all the same.
const OBSTACLE =
    LazySets.VPolygon([SVector(-0.39, 0.0), SVector(-0.35, 0.02), SVector(-0.31, 0.0)])
const X_BOX = LazySets.Hyperrectangle(;
    low = SVector(-0.6, -0.6, -0.6, -0.6),
    high = SVector(0.6, 0.6, 0.6, 0.6),
)

@testset "obstacle carving is sound (adversarial sampling)" begin
    removed = Set(RP.infeasible_cells(X_BOX, DISC.state_grid, OBSTACLE, GEOM))
    @test !isempty(removed)
    # The domain only holds the INNER cells of the box (boundary cells are not
    # states at all), so soundness is claimed — and tested — on those.
    domain_cells = Set(MP.get_pos_from_set(DISC.state_grid, X_BOX, MP.INNER))

    n_checked = 0
    n_violations = 0
    for _ in 1:20000
        θ = SVector{4}(0.6 .* (2 .* rand(RNG, 4) .- 1))
        pos = MP.get_pos_by_coord(DISC.state_grid, θ)
        (pos in domain_cells && pos ∉ removed) || continue
        n_checked += 1
        # A configuration in a kept cell never puts the swing foot inside the
        # obstacle (that is exactly what OUTER carving certifies).
        RP.swing_foot_position(θ, GEOM, true) ∈ OBSTACLE && (n_violations += 1)
    end
    @test n_checked > 0
    @test n_violations == 0
end

@testset "target cells contain the foothold manifold" begin
    foothold = SVector(0.2, 0.0)
    cells = RP.target_cells(DISC.state_grid, foothold, GEOM, X_BOX)
    @test !isempty(cells)
    # Deduplicated: positions are unique by construction of the Set.
    @test length(cells) == length(Set(cells))

    dev = RP.swing_foot_deviation_bound(GEOM, MP.get_h(DISC.state_grid) ./ 2)
    for pos in cells
        center = SVector{4}(MP.get_coord_by_pos(DISC.state_grid, pos))
        foot = RP.swing_foot_position(center, GEOM, true)
        # The cell contains a configuration with the foot exactly on the
        # foothold, so the center's foot is within the cell deviation bound.
        @test sqrt(sum((foot - foothold) .^ 2)) <= dev + 1e-9
    end
end

@testset "one footstep, end to end: exact lattice, no obstacle crossing" begin
    domain = RP.RobotDomainConfig(;
        x_lb = SVector{4}(LazySets.low(X_BOX)),
        x_ub = SVector{4}(LazySets.high(X_BOX)),
        u_lb = SVector(-DISC.u_max, -DISC.u_max, -DISC.u_max, -DISC.u_max),
        u_ub = SVector(DISC.u_max, DISC.u_max, DISC.u_max, DISC.u_max),
    )
    sys = RP.system(;
        tstep = DISC.tstep,
        domain = domain,
        obstacle = OBSTACLE,
        state_grid = DISC.state_grid,
        geometry = GEOM,
    )
    @test sys.X isa UT.SetMinus

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(sys, nothing),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), DISC.state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), DISC.input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    @test SY.get_n_input(abstract_system) == 3^4

    # Exact lattice: every transition moves by exactly one cell per axis at
    # most — the concrete displacement tstep * u is a grid vector. Violations
    # are counted (not @test-ed) to keep the record count independent of the
    # abstraction size.
    h = MP.get_h(DISC.state_grid)
    n_bad_displacement = 0
    n_bad_jump = 0
    for (target, source, symbol) in SY.enum_transitions(abstract_system)
        xs = SY.get_concrete_state(abstract_system, source)
        xt = SY.get_concrete_state(abstract_system, target)
        u = SY.get_concrete_input(abstract_system, symbol)
        maximum(abs.(xt - xs - DISC.tstep .* u)) <= 1e-9 || (n_bad_displacement += 1)
        maximum(abs.(xt - xs) ./ h) <= 1 + 1e-9 || (n_bad_jump += 1)
    end
    @test n_bad_displacement == 0
    @test n_bad_jump == 0

    x0 = SVector(0.2, 0.0, -0.2, 0.0)
    foothold = SVector(0.2, 0.0)
    step_pb = RP.step_problem(sys, DISC.state_grid; x0 = x0, foothold = foothold)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

    reached(x) = x ∈ step_pb.target_set
    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        controller,
        x0,
        100;
        stopping = reached,
    )
    xs = collect(ST.states(traj))
    @test reached(xs[end])

    removed = Set(RP.infeasible_cells(X_BOX, DISC.state_grid, OBSTACLE, GEOM))
    for k in 1:(length(xs) - 1)
        # Bisimulation witness: each concrete step lands exactly on a grid
        # point, one cell away at most.
        @test maximum(abs.(xs[k + 1] - xs[k]) ./ h) <= 1 + 1e-9
        # Inter-sample safety: the straight joint-space segment between two
        # samples never puts the swing foot in the obstacle or underground.
        # (12 points so no sample falls exactly on the λ = 0.5 cell boundary.)
        for λ in range(0.0, 1.0; length = 12)
            θ = SVector{4}((1 - λ) .* xs[k] + λ .* xs[k + 1])
            @test MP.get_pos_by_coord(DISC.state_grid, θ) ∉ removed
            foot = RP.swing_foot_position(θ, GEOM, true)
            @test foot ∉ OBSTACLE
        end
    end
end

end # module
