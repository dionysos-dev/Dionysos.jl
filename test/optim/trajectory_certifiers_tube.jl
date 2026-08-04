module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets

# Pure tube-construction helpers of the uniform-grid trajectory certifier (no abstraction
# build, so these run in the fast set). The end-to-end certifier is covered separately in
# the `:slow` `trajectory_certifiers.jl`.
const UGC = AB.UniformGridTrajectoryCertifier

@testset "uniform-grid certifier: fixed densification" begin
    xs = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0)]

    # n_between <= 0 is a no-op returning the same vector.
    @test UGC.densify_points_fixed(xs, 0) === xs

    # Insert 2 points between each consecutive pair, keeping the endpoints.
    dz = UGC.densify_points_fixed(xs, 2)
    @test length(dz) == length(xs) + (length(xs) - 1) * 2   # 3 + 2*2 = 7
    @test dz[1] == xs[1]
    @test dz[end] == xs[end]
    @test dz[2] ≈ SVector(1 / 3, 0.0)                        # evenly spaced on the segment
    @test dz[3] ≈ SVector(2 / 3, 0.0)
end

@testset "uniform-grid certifier: max-step densification" begin
    ys = [SVector(0.0, 0.0), SVector(0.05, 0.0), SVector(1.0, 0.0)]

    # Non-positive max_step is a no-op.
    @test UGC.densify_points_maxstep(ys, 0.0) === ys

    # First segment (‖Δ‖∞ = 0.05 ≤ 0.1) needs no split; the long second segment does.
    dm = UGC.densify_points_maxstep(ys, 0.1)
    @test dm[1] == ys[1]
    @test dm[end] == ys[end]
    @test ys[2] in dm                    # the short-segment endpoint is preserved
    @test length(dm) > length(ys)        # the long segment was subdivided
    @test all(maximum(abs.(dm[k + 1] - dm[k])) <= 0.1 + 1e-9 for k in 1:(length(dm) - 1))
end

@testset "uniform-grid certifier: build_tube" begin
    traj = ST.Trajectory([SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0)])

    # Fixed densification path inside build_tube (n_between > 0).
    tube = UGC.build_tube(traj, 0.2; n_between = 1, enforce_safe_max_step = true)
    @test !isnothing(tube)

    # max_step above 0.5 * radius emits a gap warning (only when the safe cap is off).
    tube_warn = @test_logs (:warn,) UGC.build_tube(
        traj,
        0.2;
        max_step = 0.5,
        enforce_safe_max_step = false,
    )
    @test !isnothing(tube_warn)

    # Intersecting with a bounded domain keeps the tube within it.
    dom = UT.box(SVector(-1.0, -1.0), SVector(2.0, 2.0))
    tube_dom = UGC.build_tube(traj, 0.2; enforce_safe_max_step = true, X_domain = dom)
    @test !isnothing(tube_dom)
end

end # module TestMain
