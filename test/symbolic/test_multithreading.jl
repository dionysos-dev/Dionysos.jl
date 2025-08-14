module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
using LinearAlgebra

sleep(0.1) # used for good printing
println("Started multithreading test")

# Test helper function to build example system
function build_test_system(; n_per_dim = 20, tstep = 1.0, input_step = 1.0)
    # 3D state domain
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)
    Xgrid = DO.GridFree(lb, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

    # 3D input domain
    lb_u = SVector(-1.0, -1.0, -1.0)
    ub_u = SVector(1.0, 1.0, 1.0)
    h_u = SVector(input_step, input_step, input_step)
    Ugrid = DO.GridFree(lb_u, h_u)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb_u, ub_u), DO.OUTER)

    # Linear dynamics: dx/dt = A*x + B*u
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
    function F_sys(x, u)
        return A * x + B * u
    end

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        3,
        3,
        nothing,
        nothing,
    )

    return Xfull, Ufull, concrete_system
end

# Test function to verify consistency between serial and threaded execution
function test_multithreading_consistency(method_name, concrete_system, Xfull, Ufull)
    # Serial execution
    sym_serial = SY.NewSymbolicModelListList(Xfull, Ufull)
    SY.compute_abstract_system_from_concrete_system!(
        sym_serial,
        concrete_system;
        verbose = false,
        threaded = false,
    )
    transitions_serial = Set(SY.enum_transitions(sym_serial.autom))
    n_transitions_serial = length(transitions_serial)

    # Threaded execution (if available)
    if Threads.nthreads() > 1
        sym_threaded = SY.NewSymbolicModelListList(Xfull, Ufull)
        SY.compute_abstract_system_from_concrete_system!(
            sym_threaded,
            concrete_system;
            verbose = false,
            threaded = true,
        )
        transitions_threaded = Set(SY.enum_transitions(sym_threaded.autom))
        n_transitions_threaded = length(transitions_threaded)

        # Verify consistency
        is_equal = transitions_serial == transitions_threaded
        return (
            consistent = is_equal,
            n_serial = n_transitions_serial,
            n_threaded = n_transitions_threaded,
        )
    else
        return (
            consistent = true,
            n_serial = n_transitions_serial,
            n_threaded = n_transitions_serial,
        )
    end
end

# Test function to measure speedup for a specific method
function measure_speedup(method_name, concrete_system, Xfull, Ufull; repeats = 3)
    if Threads.nthreads() == 1
        return 1.0  # No speedup possible with single thread
    end

    serial_times = Float64[]
    threaded_times = Float64[]

    for _ in 1:repeats
        # Serial measurement
        sym_serial = SY.NewSymbolicModelListList(Xfull, Ufull)
        GC.gc()
        t_serial = @elapsed SY.compute_abstract_system_from_concrete_system!(
            sym_serial,
            concrete_system;
            verbose = false,
            threaded = false,
        )
        push!(serial_times, t_serial)

        # Threaded measurement
        sym_threaded = SY.NewSymbolicModelListList(Xfull, Ufull)
        GC.gc()
        t_threaded = @elapsed SY.compute_abstract_system_from_concrete_system!(
            sym_threaded,
            concrete_system;
            verbose = false,
            threaded = true,
        )
        push!(threaded_times, t_threaded)
    end

    avg_serial = sum(serial_times) / length(serial_times)
    avg_threaded = sum(threaded_times) / length(threaded_times)

    return avg_serial / avg_threaded
end

@testset "Multithreading Consistency" begin
    @testset "DiscreteTimeCenteredSimulation" begin
        Xfull, Ufull, concrete_system = build_test_system(; n_per_dim = 15)

        # Create discrete system
        cont_center = ST.ContinuousTimeCenteredSimulation(concrete_system)
        discrete_system = ST.discretize(cont_center, 1.0)

        result = test_multithreading_consistency(
            "CenteredSimulation",
            discrete_system,
            Xfull,
            Ufull,
        )

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        # Measure and display speedup
        speedup = measure_speedup("CenteredSimulation", discrete_system, Xfull, Ufull)
        println("CenteredSimulation: $(round(speedup, digits=2))× speedup")
    end

    @testset "DiscreteTimeSystemOverApproximation" begin
        Xfull, Ufull, concrete_system = build_test_system(; n_per_dim = 15)

        # Create over-approximation system
        function simple_over_approx(elem, u)
            center = UT.get_center(elem)
            radius = UT.get_r(elem)
            new_radius = radius * 1.1 .+ 0.01
            F_sys = concrete_system.f
            new_center = F_sys(center, u)
            return UT.HyperRectangle(new_center - new_radius, new_center + new_radius)
        end
        discrete_system = ST.discretize_continuous_system(concrete_system, 1.0)
        over_approx_system =
            ST.DiscreteTimeOverApproximationMap(discrete_system, simple_over_approx)

        result = test_multithreading_consistency(
            "OverApproximation",
            over_approx_system,
            Xfull,
            Ufull,
        )

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        # Measure and display speedup
        speedup = measure_speedup("OverApproximation", over_approx_system, Xfull, Ufull)
        println("OverApproximation: $(round(speedup, digits=2))× speedup")
    end

    @testset "DiscreteTimeGrowthBound" begin
        Xfull, Ufull, concrete_system = build_test_system(; n_per_dim = 15)

        # Create growth bound system
        function growth_bound_map(r, u)
            return r * (1.0 + 0.1 * norm(u)) .+ 0.01
        end
        discrete_system = ST.discretize_continuous_system(concrete_system, 1.0)
        growth_bound_system = ST.DiscreteTimeGrowthBound(discrete_system, growth_bound_map)

        result = test_multithreading_consistency(
            "GrowthBound",
            growth_bound_system,
            Xfull,
            Ufull,
        )

        @test result.consistent == true
        @test result.n_serial > 0
        if Threads.nthreads() > 1
            @test result.n_serial == result.n_threaded
        end

        # Measure and display speedup
        speedup = measure_speedup("GrowthBound", growth_bound_system, Xfull, Ufull)
        println("GrowthBound: $(round(speedup, digits=2))× speedup")
    end
end

sleep(0.1) # used for good printing
println("End Test")

end # module TestMain
