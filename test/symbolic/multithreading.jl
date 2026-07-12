module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
using LinearAlgebra

# ----------------------------
# Build finite mappings + system (Mapping-based)
# ----------------------------
function build_test_system(;
    n_per_dim::Int = 20,
    tstep::Float64 = 1.0,
    input_step::Float64 = 1.0,
)
    # ---- X mapping: [0,1]^3 sampled with n_per_dim points each axis ----
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)

    Xgrid = MP.GridFree(lb, h)
    Xmap = MP.ExplicitGridMapping(Xgrid)

    # positions 0:(n_per_dim-1) in each dim
    for i in 0:(n_per_dim - 1), j in 0:(n_per_dim - 1), k in 0:(n_per_dim - 1)
        MP.add_pos!(Xmap, (i, j, k))
    end

    # ---- U mapping: [-1,1]^3 with step input_step ----
    lb_u = SVector(-1.0, -1.0, -1.0)
    ub_u = SVector(1.0, 1.0, 1.0)
    h_u = SVector(input_step, input_step, input_step)

    Ugrid = MP.GridFree(lb_u, h_u)
    Umap = MP.ExplicitGridMapping(Ugrid)

    # positions that hit coords in [lb_u, ub_u] on this grid
    # coord = lb_u + pos*h_u, so pos ranges 0..round((ub-lb)/h)
    nu = ntuple(d -> Int(round((ub_u[d] - lb_u[d]) / h_u[d])), 3)
    for i in 0:nu[1], j in 0:nu[2], k in 0:nu[3]
        MP.add_pos!(Umap, (i, j, k))
    end

    # ---- Linear dynamics: dx/dt = A*x + B*u ----
    A = @SMatrix [
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        -1.0 -1.0 -1.0
    ]
    B = @SMatrix [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0
    ]
    F_sys(x, u) = A * x + B * u

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        3,
        3,
        nothing,
        nothing,
    )

    return Xmap, Umap, concrete_system
end

# ----------------------------
# Consistency check: serial vs threaded
# ----------------------------
function test_multithreading_consistency(method_name, concrete_system, Xmap, Umap)
    # Serial
    sym_serial = SY.SymbolicModelList(Xmap, Umap)
    SY.compute_abstract_system_from_concrete_system!(
        sym_serial,
        concrete_system;
        threaded = false,
    )
    transitions_serial = Set(SY.enum_transitions(sym_serial.autom))
    n_serial = length(transitions_serial)

    # Threaded (if available)
    if Threads.nthreads() > 1
        sym_threaded = SY.SymbolicModelList(Xmap, Umap)
        SY.compute_abstract_system_from_concrete_system!(
            sym_threaded,
            concrete_system;
            threaded = true,
        )
        transitions_threaded = Set(SY.enum_transitions(sym_threaded.autom))
        n_threaded = length(transitions_threaded)

        return (
            consistent = (transitions_serial == transitions_threaded),
            n_serial = n_serial,
            n_threaded = n_threaded,
        )
    else
        return (consistent = true, n_serial = n_serial, n_threaded = n_serial)
    end
end

# ----------------------------
# Speed measurement
# ----------------------------
function measure_speedup(method_name, concrete_system, Xmap, Umap; repeats::Int = 3)
    Threads.nthreads() == 1 && return 1.0

    serial_times = Float64[]
    threaded_times = Float64[]

    for _ in 1:repeats
        sym_serial = SY.SymbolicModelList(Xmap, Umap)
        GC.gc()
        t_serial = @elapsed SY.compute_abstract_system_from_concrete_system!(
            sym_serial,
            concrete_system;
            threaded = false,
        )
        push!(serial_times, t_serial)

        sym_threaded = SY.SymbolicModelList(Xmap, Umap)
        GC.gc()
        t_threaded = @elapsed SY.compute_abstract_system_from_concrete_system!(
            sym_threaded,
            concrete_system;
            threaded = true,
        )
        push!(threaded_times, t_threaded)
    end

    avg_serial = sum(serial_times) / length(serial_times)
    avg_threaded = sum(threaded_times) / length(threaded_times)

    return avg_serial / avg_threaded
end

# ----------------------------
# Tests
# ----------------------------
@testset "Multithreading Consistency" begin
    @testset "DiscreteTimeCenteredSimulation" begin
        Xmap, Umap, concrete_system = build_test_system(; n_per_dim = 15)

        cont_center = ST.ContinuousTimeCenteredSimulation(concrete_system)
        discrete_system = ST.discretize(cont_center, 1.0)

        result = test_multithreading_consistency(
            "CenteredSimulation",
            discrete_system,
            Xmap,
            Umap,
        )

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        speedup = measure_speedup("CenteredSimulation", discrete_system, Xmap, Umap)
        println("CenteredSimulation: $(round(speedup, digits=2))× speedup")
    end

    @testset "DiscreteTimeSystemOverApproximation" begin
        Xmap, Umap, concrete_system = build_test_system(; n_per_dim = 15)

        function simple_over_approx(elem, u)
            center = LazySets.center(elem)
            radius = LazySets.radius_hyperrectangle(elem)
            new_radius = radius * 1.1 .+ 0.01
            F_sys = concrete_system.f
            new_center = F_sys(center, u)
            return UT.box(new_center - new_radius, new_center + new_radius)
        end

        discrete_system = ST.discretize_continuous_system(concrete_system, 1.0)
        over_approx_system =
            ST.DiscreteTimeOverApproximationMap(discrete_system, simple_over_approx)

        result = test_multithreading_consistency(
            "OverApproximation",
            over_approx_system,
            Xmap,
            Umap,
        )

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        speedup = measure_speedup("OverApproximation", over_approx_system, Xmap, Umap)
        println("OverApproximation: $(round(speedup, digits=2))× speedup")
    end

    @testset "DiscreteTimeGrowthBound" begin
        Xmap, Umap, concrete_system = build_test_system(; n_per_dim = 15)

        growth_bound_map(r, u) = r * (1.0 + 0.1 * norm(u)) .+ 0.01

        discrete_system = ST.discretize_continuous_system(concrete_system, 1.0)
        growth_bound_system = ST.DiscreteTimeGrowthBound(discrete_system, growth_bound_map)

        result =
            test_multithreading_consistency("GrowthBound", growth_bound_system, Xmap, Umap)

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        speedup = measure_speedup("GrowthBound", growth_bound_system, Xmap, Umap)
        println("GrowthBound: $(round(speedup, digits=2))× speedup")
    end
end

println("End Test")

end # module TestMain
