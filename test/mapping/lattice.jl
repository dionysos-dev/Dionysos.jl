module TestLattice

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI
import Logging

@testset "intersample_safe_time_step" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.2))

    @test MP.intersample_safe_time_step(grid, SVector(1.0, 1.0)) ≈ 0.1
    @test MP.intersample_safe_time_step(grid, SVector(2.0, 0.5)) ≈ 0.05
    # A zero velocity bound leaves that axis unconstrained.
    @test MP.intersample_safe_time_step(grid, SVector(0.0, 0.5)) ≈ 0.4
    @test_throws ErrorException MP.intersample_safe_time_step(grid, SVector(1.0))
end

@testset "is_lattice_exact" begin
    state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1))

    # tstep * du = 0.1 = hx: exact.
    exact_input = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    @test MP.is_lattice_exact(state_grid, exact_input, 0.1)
    # Multiples of the cell width are exact too.
    @test MP.is_lattice_exact(state_grid, exact_input, 0.2)

    # tstep * du = 0.07: not on the lattice.
    @test !MP.is_lattice_exact(state_grid, exact_input, 0.07)

    # Input-grid origin off the lattice breaks exactness.
    offset_input = MP.GridFree(SVector(0.3, 0.0), SVector(1.0, 1.0))
    @test !MP.is_lattice_exact(state_grid, offset_input, 0.1)

    # Dimension mismatch is never exact.
    @test !MP.is_lattice_exact(state_grid, MP.GridFree(SVector(0.0), SVector(1.0)), 0.1)
end

# The end-to-end warning: a 1-D integrator over a domain with a hole thinner
# than one step's displacement. With tstep = 0.3 and |u| ≤ 1 on a 0.1 grid, a
# step moves up to 3 cells: the abstraction jumps the hole and must say so.
function _integrator_optimizer(; tstep)
    sys = single_integrator(; n = 1, xbound = 1.0, ubound = 1.0)
    hole = LazySets.Hyperrectangle(; low = SVector(-0.549), high = SVector(-0.351))
    X = UT.set_minus(sys.X, hole)
    carved = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        sys.f,
        1,
        1,
        X,
        sys.U,
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(carved, X),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0), SVector(0.1)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(1.0)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), u -> SMatrix{1, 1}(0.0))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    return optimizer
end

@testset "inter-sample jump warning" begin
    optimizer = _integrator_optimizer(; tstep = 0.3)
    @test_logs (:warn, r"jump.*cells per step") match_mode = :any MOI.optimize!(optimizer)
end

@testset "intersample_checked silences the warning" begin
    optimizer = _integrator_optimizer(; tstep = 0.3)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("intersample_checked"), true)
    @test_logs min_level = Logging.Warn MOI.optimize!(optimizer)
end

@testset "cells_on_segment" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1))

    # Degenerate segment: the single containing cell.
    @test MP.cells_on_segment(grid, SVector(0.0, 0.0), SVector(0.0, 0.0)) == [(0, 0)]

    # Axis moves of one and two cells (cell centers to cell centers).
    @test MP.cells_on_segment(grid, SVector(0.0, 0.0), SVector(0.1, 0.0)) ==
          [(0, 0), (1, 0)]
    @test MP.cells_on_segment(grid, SVector(0.0, 0.0), SVector(0.2, 0.0)) ==
          [(0, 0), (1, 0), (2, 0)]

    # A (2, 1)-cell move visits the staircase cells in order.
    @test MP.cells_on_segment(grid, SVector(0.0, 0.0), SVector(0.2, 0.1)) ==
          [(0, 0), (1, 0), (1, 1), (2, 1)]

    # An exact corner crossing yields the two diagonal cells (the corner point
    # lies in the closure of both).
    @test MP.cells_on_segment(grid, SVector(0.0, 0.0), SVector(0.1, 0.1)) ==
          [(0, 0), (1, 1)]

    # Every sampled point of a random segment lies in the closure of a
    # returned cell (soundness), for segments not starting on cell centers.
    a, b = SVector(0.031, -0.017), SVector(-0.24, 0.113)
    cells = MP.cells_on_segment(grid, a, b)
    for λ in range(0.0, 1.0; length = 101)
        p = (1 - λ) .* a + λ .* b
        @test any(all(abs(p[i] - 0.1 * c[i]) <= 0.05 + 1e-12 for i in 1:2) for c in cells)
    end
end

end # module
