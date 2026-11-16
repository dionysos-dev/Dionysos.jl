# Tests for SampledAbstraction.jl, on the simple pendulum (2 states) and the double pendulum (4 states)
# of `problems/Pendulum`, so that every abstraction takes a fraction of a second. Run with an
# environment that has Dionysos, Symbolics, IntervalArithmetic, CoreMath, MathematicalSystems, LazySets,
# StaticArrays, JLD2, MathOptInterface and Plots (the pendulum problems load it):
#   julia --project=<env> tests/test_sampled_abstraction.jl
using Test
using Random
using StaticArrays
using Symbolics
using IntervalArithmetic          # loads Dionysos' interval/Symbolics extension
import CoreMath
import Dionysos
import LazySets
import MathematicalSystems as MS
import MathOptInterface as MOI
import JLD2

const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MP = Dionysos.Mapping
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction
const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")

include(joinpath(PROBLEMS, "Pendulum", "simple_pendulum.jl"))
include(joinpath(PROBLEMS, "Pendulum", "double_pendulum.jl"))
include(joinpath(@__DIR__, "..", "LocalGrowthBound.jl"))
include(joinpath(@__DIR__, "..", "SampledAbstraction.jl"))   # includes BipedAbstraction.jl
using .LocalGrowthBound: SymbolicDynamics, cellwise_bound
using .SampledAbstraction
import .SampledAbstraction: _mix64, _stream, _unit, _cell_key, _uniform_design

const BA = BipedAbstraction

# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------

abstract_system(opt) = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
transitions(opt) = Set(SY.enum_transitions(abstract_system(opt)))

# The targets of every enabled (source, input) pair.
function successors(opt)
    d = Dict{Tuple{Int, Int}, Set{Int}}()
    for (target, source, input) in SY.enum_transitions(abstract_system(opt))
        push!(get!(d, (source, input), Set{Int}()), target)
    end
    return d
end

# `inner` is an abstraction inside `outer`: wherever `outer` enables a pair, `inner` does too, and its
# targets are among those of `outer`.
function inside(inner, outer)
    for (pair, targets) in outer
        haskey(inner, pair) || return false
        issubset(inner[pair], targets) || return false
    end
    return true
end

# A sampler that counts its calls, to see whether a build computes anything. Its `string` is
# constant: that is what identifies a sampler in the provenance.
struct CountingSampler
    calls::Base.RefValue{Int}
end
(s::CountingSampler)(cell, n, key) = (s.calls[] += 1; RandomSampler()(cell, n, key))
Base.show(io::IO, ::CountingSampler) = print(io, "CountingSampler()")

# ---------------------------------------------------------------------------------------------
# The pendulum
# ---------------------------------------------------------------------------------------------

pendulum = SimplePendulum.system()
tstep = 0.1
state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.5))
input_grid = MP.GridFree(SVector(0.0), SVector(3.0))
grids = (; state_grid, input_grid)

sampled(system = pendulum; kwargs...) =
    build_sampled_abstraction(system; grids..., time_step = tstep, kwargs...)

@testset "SampledAbstraction" begin
    Random.seed!(1234)

    @testset "the generator is SplitMix64, the same everywhere" begin
        # The reference outputs of SplitMix64 started from 0 (Steele, Lea, Flood; Vigna's C code).
        @test [_stream(UInt64(0), i) for i in 0:2] == [0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4, 0x06c45d188009454f]
        u = [_unit(UInt64(7), i) for i in 0:9999]
        @test all(0 .<= u .< 1)
        @test _unit(UInt64(7), 5) == u[6]                    # a pure function of the position
        @test abs(sum(u) / length(u) - 0.5) < 5 / sqrt(12 * length(u))
        # Seeds and cells get their own streams.
        c = SVector(0.3, -0.5)
        @test _cell_key(0, c) == _cell_key(0, c)
        @test _cell_key(0, c) != _cell_key(1, c)
        @test _cell_key(0, c) != _cell_key(0, SVector(0.3, -0.0))
        @test _cell_key(0, SVector(0.0, 0.0)) == _cell_key(0, SVector(-0.0, 0.0))
    end

    @testset "random sampler" begin
        cell = LazySets.Hyperrectangle(;
            low = SVector(1.0, -2.0, 0.0),
            high = SVector(1.5, 2.0, 3.0),
        )
        c, ρ = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        sampler = RandomSampler()
        key = _cell_key(3, SVector{3, Float64}(c))
        pts = sampler(cell, 20000, key)
        @test length(pts) == 20000 && eltype(pts) == SVector{3, Float64}
        @test all(p -> p in cell, pts)
        # Reproducible, and a prefix does not depend on how many points are asked for.
        @test pts == sampler(cell, 20000, key)
        @test sampler(cell, 7, key) == pts[1:7]
        @test sampler(cell, 7, key + 1) != pts[1:7]
        # Independent and uniform: mean, variance, halves and correlations, within 5 sigma.
        N = length(pts)
        for j in 1:3
            x = [(p[j] - c[j]) / ρ[j] for p in pts]               # uniform in [-1, 1)
            @test abs(sum(x) / N) < 5 * sqrt(1 / 3 / N)
            @test isapprox(sum(abs2, x) / N, 1 / 3; rtol = 0.05)
            @test abs(count(>(0), x) - N / 2) < 5 * sqrt(N / 4)
            for k in (j + 1):3
                y = [(p[k] - c[k]) / ρ[k] for p in pts]
                @test abs(sum(x .* y) / N / (1 / 3)) < 5 / sqrt(N)
            end
        end
        # A cell of a Dionysos grid has SVector storage, but any hyperrectangle will do.
        generic = LazySets.Hyperrectangle([1.25, 0.0, 1.5], [0.25, 2.0, 1.5])
        @test all(p -> p in generic, sampler(generic, 50, key))
    end

    @testset "the uniform design: center, corners, finer lattices" begin
        points(M) = Set(Tuple(c) for c in eachcol(M))
        # The lattice with `k` points per side of the cube [-1, 1]ⁿ (k - 1 is a power of 2: exact).
        lattice(n, k) = Set(
            Iterators.product(
                ntuple(_ -> [-1.0 + 2.0 * i / (k - 1) for i in 0:(k - 1)], n)...,
            ),
        )
        corners(n) = Set(Iterators.product(ntuple(_ -> (-1.0, 1.0), n)...))

        # One point is the center; the center and the 2ⁿ corners are 17 points in four dimensions.
        @test _uniform_design(4, 1) == zeros(4, 1)
        design = _uniform_design(4, 17)
        @test size(design) == (4, 17) && iszero(design[:, 1])
        @test points(design[:, 2:end]) == corners(4)
        @test points(_uniform_design(2, 5)) == union(corners(2), Set([(0.0, 0.0)]))
        @test points(_uniform_design(1, 3)) == Set([(0.0,), (-1.0,), (1.0,)])
        # Every level that is complete is a regular lattice: 3 points per side, then 5, 9.
        for (n, k) in
            ((1, 3), (1, 5), (1, 9), (2, 3), (2, 5), (2, 9), (3, 3), (3, 5), (4, 3), (4, 5))
            @test points(_uniform_design(n, k^n)) == lattice(n, k)
        end
        # Anything else is the beginning of the next lattice: distinct points of the cube, and the
        # design of `N` points is the first `N` of that of more (the abstraction only grows with N).
        D = _uniform_design(3, 100)
        @test length(points(D)) == 100 && all(abs.(D) .<= 1)
        @test all(N -> _uniform_design(3, N) == D[:, 1:N], (1, 2, 8, 9, 10, 27, 28, 50, 99))
        # A partial level is spread over the cube. First the corners come in opposite pairs, so
        # that the center and two corners are on a diagonal ...
        two = _uniform_design(2, 3)
        @test two[:, 3] == -two[:, 2]
        @test _uniform_design(4, 3)[:, 3] == -_uniform_design(4, 3)[:, 2]
        # ... and so do the midpoints of a square, or of an interval.
        @test _uniform_design(2, 7)[:, 7] == -_uniform_design(2, 7)[:, 6]
        @test _uniform_design(1, 7)[:, 6:7] == [-0.75 0.75]
        # Each point is a farthest one: of the points of its level that are not chosen yet, none is
        # farther from the points selected before it (checked by brute force).
        function farthest_first(n, N, level_points, first)
            D = _uniform_design(n, N)
            ok = true
            for t in first:N
                dmin(q) = minimum(sum(abs2, q .- D[:, s]) for s in 1:(t - 1))
                taken = Set(Tuple(D[:, s]) for s in 1:(t - 1))
                rest = [q for q in level_points if !(q in taken)]
                ok &= dmin(collect(D[:, t])) >= maximum(dmin(collect(q)) for q in rest)
            end
            return ok
        end
        mid(n) = setdiff(lattice(n, 3), corners(n), Set([ntuple(_ -> 0.0, n)]))
        @test farthest_first(2, 9, mid(2), 6)
        @test farthest_first(3, 27, mid(3), 10)
        @test farthest_first(2, 25, setdiff(lattice(2, 5), lattice(2, 3)), 10)
        # A design that would need an enormous lattice is refused, and says what to use instead.
        err = try
            _uniform_design(10, 100_000)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin(":random", err.msg)
        @test_throws ArgumentError _uniform_design(2, 0)

        # The sampler puts it in a cell: the center, then the corners.
        cell = LazySets.Hyperrectangle(;
            low = SVector(1.0, -2.0, 0.0, 5.0),
            high = SVector(1.5, 2.0, 3.0, 6.0),
        )
        c, ρ = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        pts = UniformSampler()(cell, 17, UInt64(0))
        @test pts[1] == c && length(pts) == 17
        @test Set(Tuple.(pts[2:end])) ==
              Set(Iterators.product(ntuple(j -> (c[j] - ρ[j], c[j] + ρ[j]), 4)...))
        @test all(p -> p in cell, pts)
        # No key, no seed. Prepared once for a dimension and a number of points, or not.
        prepared = SampledAbstraction.prepare(UniformSampler(), 4, 17)
        @test prepared(cell, 17, UInt64(12345)) == pts
        @test_throws ArgumentError prepared(cell, 18, UInt64(0))
        @test string(prepared) == string(UniformSampler()) == "UniformSampler()"
        @test string(RandomSampler()) == "RandomSampler()"
    end

    @testset "the image of a cell is looked up with the quantizer of the grid" begin
        opt = sampled(; nsamples = 4)
        sym = abstract_system(opt)
        mapping = SY.get_state_mapping(sym)
        lookup(points) = MP.get_states_from_set_strict(
            mapping,
            SampledImage([SVector{2, Float64}(p) for p in points]),
            MP.OUTER,
        )
        a, b = SVector(0.3, 0.5), SVector(-0.9, 1.0)
        qa, qb = SY.get_abstract_state(sym, a), SY.get_abstract_state(sym, b)
        @test qa !== nothing && qb !== nothing && qa != qb
        @test lookup([a, b, a]) == ([qa, qb], true)                 # the union, each cell once
        @test lookup([a, SVector(0.31, 0.52)]) == ([qa], true)      # two points of a cell
        # A point on the edge of two cells goes where the quantizer puts it (not to both).
        edge = SVector(0.15, 0.0)
        @test lookup([edge]) == ([SY.get_abstract_state(sym, edge)], true)
        # A point outside the domain, or that is not a number, disables the pair.
        for bad in
            (SVector(10.0, 0.0), SVector(NaN, 0.0), SVector(0.0, Inf), SVector(1e300, 0.0))
            targets, allin = lookup([a, bad])
            @test !allin
        end
        img = SampledImage([a, b])
        @test LazySets.dim(img) == 2
        @test LazySets.low(img) == SVector(-0.9, 0.5) &&
              LazySets.high(img) == SVector(0.3, 1.0)
        # The same with an implicit mapping.
        implicit = abstract_system(sampled(; nsamples = 4, use_implicit_mapping = true))
        imapping = SY.get_state_mapping(implicit)
        @test MP.get_states_from_set_strict(imapping, SampledImage([a, b]), MP.OUTER) == (
            [SY.get_abstract_state(implicit, a), SY.get_abstract_state(implicit, b)],
            true,
        )
        @test !last(
            MP.get_states_from_set_strict(
                imapping,
                SampledImage([SVector(10.0, 0.0)]),
                MP.OUTER,
            ),
        )
    end

    @testset "the abstraction is reproducible" begin
        a = sampled(; nsamples = 16, seed = 1)
        @test !isempty(transitions(a))
        @test transitions(a) == transitions(sampled(; nsamples = 16, seed = 1))
        @test transitions(a) != transitions(sampled(; nsamples = 16, seed = 2))
        # ... whatever the order of the cells and the threads.
        threaded =
            sampled(; nsamples = 16, seed = 1, execution_backend = SY.ThreadedBackend(0.2))
        @test transitions(threaded) == transitions(a)
        # It needs only the samples of the cell: `N` samples are the first `N` of more.
        small, large =
            successors(sampled(; nsamples = 4)), successors(sampled(; nsamples = 32))
        @test inside(small, large)
        @test length(transitions(sampled(; nsamples = 4))) <
              length(transitions(sampled(; nsamples = 32)))
    end

    @testset "uniform sampling" begin
        uniform(N; kwargs...) = sampled(; nsamples = N, sampling = :uniform, kwargs...)
        center = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.CENTER_SIMULATION,
        )
        # One point is the center simulation, and `sampling` may be written as a string.
        @test transitions(uniform(1)) == transitions(center)
        @test transitions(sampled(; nsamples = 9, sampling = "uniform")) ==
              transitions(uniform(9))

        # The center and the four corners of every cell, simulated by hand with the RK4 map of
        # Dionysos: the cells of the images, and nothing if one of them leaves the domain.
        opt = uniform(5)
        sym = abstract_system(opt)
        found = successors(opt)
        F = ST.get_system_map(SampledSimulation(pendulum, tstep; nsamples = 1))
        mismatches = 0
        for s in 1:SY.get_n_state(sym), i in 1:SY.get_n_input(sym)
            cell, u = SY.get_concrete_elem(sym, s), SY.get_concrete_input(sym, i)
            c, ρ = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
            xs = [
                c,
                c .+ ρ .* SVector(1, 1),
                c .+ ρ .* SVector(1, -1),
                c .+ ρ .* SVector(-1, 1),
                c .+ ρ .* SVector(-1, -1),
            ]
            targets = [SY.get_abstract_state(sym, F(x, u)) for x in xs]
            expected =
                all(q -> SY.is_allowed_state(sym, q), targets) ? Set(targets) : nothing
            mismatches += get(found, (s, i), nothing) == expected ? 0 : 1
        end
        @test mismatches == 0
        @test length(found) > 100

        # No random number: another seed is the same abstraction, and so are the threads.
        @test transitions(uniform(9; seed = 1)) == transitions(uniform(9; seed = 2))
        @test transitions(uniform(9)) ==
              transitions(uniform(9; execution_backend = SY.ThreadedBackend(0.2)))
        # The more points, the more successors: nothing is lost on the way to a finer lattice.
        @test inside(successors(uniform(5)), successors(uniform(9)))
        @test inside(successors(uniform(9)), successors(uniform(25)))
        @test length(transitions(uniform(5))) < length(transitions(uniform(25)))
        # It is not the random sampling, which has no corners.
        @test transitions(uniform(25)) != transitions(sampled(; nsamples = 25, seed = 1))

        # What the options mean for the provenance: the seed is not part of a uniform design.
        design(; kwargs...) = SampledAbstraction.get_provenance(
            SampledSimulation(pendulum, tstep; nsamples = 9, kwargs...),
        )
        @test design(; sampling = :uniform, seed = 1) ==
              design(; sampling = :uniform, seed = 2)
        @test design(; sampling = :uniform).seed === nothing
        @test design(; sampling = :random, seed = 1) !=
              design(; sampling = :random, seed = 2)
        @test design(; sampling = :random, seed = 1) !=
              design(; sampling = :uniform, seed = 1)
        # The center is already among the points; the flag is a name, or a sampler.
        @test_throws ArgumentError SampledSimulation(
            pendulum,
            tstep;
            nsamples = 9,
            sampling = :uniform,
            include_center = true,
        )
        @test_throws ArgumentError SampledSimulation(
            pendulum,
            tstep;
            nsamples = 0,
            sampling = :uniform,
        )
        @test_throws ArgumentError SampledSimulation(
            pendulum,
            tstep;
            nsamples = 4,
            sampling = :grid,
        )
    end

    @testset "it is what Dionysos computes with the same points" begin
        # No sample and the center: the center simulation of Dionysos, transition by transition.
        center = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.CENTER_SIMULATION,
        )
        only_center = sampled(; nsamples = 0, include_center = true)
        @test !isempty(transitions(center))
        @test transitions(only_center) == transitions(center)
        # ... which is also what one sample does when it is the center: `nsamples` samples and the
        # center contain the center simulation.
        @test inside(
            successors(center),
            successors(sampled(; nsamples = 3, include_center = true)),
        )

        # The same estimator as `RANDOM_SIMULATION` (points drawn from the global generator, so
        # only the number of transitions can be compared).
        random = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.RANDOM_SIMULATION,
            n_samples = 16,
        )
        @test isapprox(
            length(transitions(sampled(; nsamples = 16))),
            length(transitions(random));
            rtol = 0.03,
        )

        # A discrete-time system is the same thing: the RK4 map of the continuous one.
        discrete = ST.discretize_continuous_system(pendulum, tstep)
        @test transitions(
            build_sampled_abstraction(discrete; grids..., nsamples = 16, seed = 1),
        ) == transitions(sampled(; nsamples = 16, seed = 1))
    end

    @testset "the successors are inside those of the over-approximations" begin
        # Every point of a cell stays in the reachable set, which the three over-approximations
        # contain: the points that are simulated cannot lead where they do not go.
        points = successors(sampled(; nsamples = 64, seed = 5, include_center = true))
        @test length(points) > 100
        global_bound = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.GROWTH,
            jacobian_bound = SimplePendulum.jacobian_bound(),
        )
        local_bound = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.GROWTH,
            jacobian_bound = cellwise_bound(SymbolicDynamics(pendulum)),
        )
        certified =
            BA.build_abstraction(BA.CertifiedRK4GrowthBound(pendulum, tstep); grids...)
        lattice = successors(sampled(; nsamples = 25, sampling = :uniform))
        for over in (global_bound, local_bound, certified)
            @test length(successors(over)) > 100              # not an empty abstraction
            @test inside(points, successors(over))
            @test inside(lattice, successors(over))           # ... nor are the corners outside it
        end
        # The bounds get tighter from the global to the local growth bound and the certified RK4.
        @test length(transitions(certified)) <=
              length(transitions(local_bound)) <=
              length(transitions(global_bound))
    end

    @testset "provenance identifies how the successors were computed" begin
        approx(; kwargs...) =
            SampledSimulation(pendulum, tstep; nsamples = 16, seed = 1, kwargs...)
        base = SampledAbstraction.get_provenance(approx())
        @test base == SampledAbstraction.get_provenance(approx())
        @test BA.provenance(approx()) == base
        @test base.approximation == "SampledSimulation" &&
              base.nsamples == 16 &&
              base.seed == 1
        other = Dict(
            "nsamples" => SampledSimulation(pendulum, tstep; nsamples = 17, seed = 1),
            "seed" => approx(; seed = 2),
            "include_center" => approx(; include_center = true),
            "sampler" => approx(; sampling = CountingSampler(Ref(0))),
            "sampling" => approx(; sampling = :uniform),
            "tstep" => SampledSimulation(pendulum, 0.2; nsamples = 16, seed = 1),
            "nsystem" => approx(; nsystem = 6),
            "dynamics" => SampledSimulation(
                SimplePendulum.system(; params = SimplePendulum.Params(; g = 9.0)),
                tstep;
                nsamples = 16,
                seed = 1,
            ),
        )
        for (what, other_approx) in other
            @test SampledAbstraction.get_provenance(other_approx) != base
        end
        # The dynamics are told apart by what they compute, not by how they are written.
        rewritten = MS.ConstrainedBlackBoxControlContinuousSystem(
            (x, u) -> SVector(x[2], u[1] - 9.81 * sin(x[1])),
            2,
            1,
            pendulum.X,
            pendulum.U,
        )
        @test SampledSimulation(rewritten, tstep; nsamples = 16, seed = 1).dynamics ==
              approx().dynamics
        # ... and a last-bit difference does not matter.
        nudged = MS.ConstrainedBlackBoxControlContinuousSystem(
            (x, u) -> SVector(x[2], u[1] - 9.81 * sin(x[1])) .* (1 + 1e-15),
            2,
            1,
            pendulum.X,
            pendulum.U,
        )
        @test SampledSimulation(nudged, tstep; nsamples = 16, seed = 1).dynamics ==
              approx().dynamics
        # What is stored in a file does not change with the session, the machine or the Julia
        # version. These values are pinned: if a change of the code moves them, every cache file and
        # chunk written before is refused (or, for the cell key, its samples are other ones), so bump
        # `SampledAbstraction.KERNEL_VERSION` with it.
        @test approx().dynamics == "fe49112e42c57eeb"
        @test _cell_key(1, SVector(0.3, 0.5)) == UInt64(16408875305039255091)
    end

    @testset "cache file and Slurm chunks" begin
        dir = mktempdir()
        file = joinpath(dir, "sampled.jld2")
        info = (; tstep, h = [0.3, 0.5])
        counter = CountingSampler(Ref(0))
        build(; kwargs...) = sampled(;
            nsamples = 16,
            seed = 1,
            sampling = counter,
            abstraction_file = file,
            abstraction_info = info,
            kwargs...,
        )
        computed = build()
        @test isfile(file)
        calls = counter.calls[]
        @test calls > 0
        loaded = build()
        @test counter.calls[] == calls                         # nothing was computed
        @test transitions(loaded) == transitions(computed)
        # The file is refused by any other construction, whatever `abstraction_info` says.
        for kwargs in (
            (; seed = 2),
            (; nsamples = 17),
            (; include_center = true),
            (; sampling = :random),
            (; sampling = :uniform),
        )
            err = try
                build(; kwargs...)
                nothing
            catch e
                e
            end
            @test err isa ErrorException && occursin("another construction", err.msg)
        end
        @test_throws ErrorException build(; time_step = 0.2)
        other_g = SimplePendulum.system(; params = SimplePendulum.Params(; g = 9.0))
        @test_throws ErrorException build_sampled_abstraction(
            other_g;
            grids...,
            time_step = tstep,
            nsamples = 16,
            seed = 1,
            sampling = counter,
            abstraction_file = file,
            abstraction_info = info,
        )
        # A uniform design has no seed: its file is the same whatever the seed, and is not the file
        # of the random sampling.
        uniform_file = joinpath(dir, "uniform.jld2")
        uniform(; kwargs...) = sampled(;
            nsamples = 9,
            sampling = :uniform,
            abstraction_file = uniform_file,
            abstraction_info = info,
            kwargs...,
        )
        written = uniform(; seed = 1)
        @test isfile(uniform_file)
        @test transitions(uniform(; seed = 2)) == transitions(written)
        @test_throws ErrorException uniform(; nsamples = 25)
        @test_throws ErrorException sampled(;
            nsamples = 9,
            seed = 1,
            abstraction_file = uniform_file,
            abstraction_info = info,
        )
        # A file of transitions computed without a record of how (Dionysos' `RANDOM_SIMULATION`) is
        # refused, and so is a sampled file by an approximation that has no record.
        plain = joinpath(dir, "plain.jld2")
        random(; kwargs...) = BA.build_abstraction(
            pendulum;
            grids...,
            time_step = tstep,
            approx_mode = UGA.RANDOM_SIMULATION,
            n_samples = 16,
            abstraction_file = plain,
            abstraction_info = info,
            kwargs...,
        )
        random()
        @test isfile(plain)
        err = try
            build_sampled_abstraction(
                pendulum;
                grids...,
                time_step = tstep,
                nsamples = 16,
                seed = 1,
                abstraction_file = plain,
                abstraction_info = info,
            )
            nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("no provenance record", err.msg)
        @test_throws ErrorException random(; abstraction_file = file)

        # The chunks of a Slurm array, merged, are the abstraction.
        full = sampled(; nsamples = 16, seed = 1)
        chunks = mktempdir()
        nchunks = 3
        chunk(id; kwargs...) = sampled(;
            nsamples = 16,
            seed = 1,
            execution_backend = SY.SlurmArrayBackend(
                nchunks,
                id,
                chunks,
                :contiguous,
                true,
            ),
            kwargs...,
        )
        foreach(chunk, 1:nchunks)
        merged =
            sampled(; nsamples = 16, seed = 1, merge_chunks = (; outdir = chunks, nchunks))
        @test transitions(merged) == transitions(full)
        # A chunk of another construction is refused: the chunks must come from one seed.
        chunk(2; seed = 2)
        @test_throws ErrorException sampled(;
            nsamples = 16,
            seed = 1,
            merge_chunks = (; outdir = chunks, nchunks),
        )
        chunk(2)
        merged =
            sampled(; nsamples = 16, seed = 1, merge_chunks = (; outdir = chunks, nchunks))
        @test transitions(merged) == transitions(full)
    end

    @testset "arguments" begin
        @test_throws ArgumentError build_sampled_abstraction(
            pendulum;
            grids...,
            nsamples = 4,
        )
        @test_throws ArgumentError build_sampled_abstraction(
            ST.discretize_continuous_system(pendulum, tstep);
            grids...,
            nsamples = 4,
            time_step = tstep,
        )
        @test_throws ArgumentError sampled(; nsamples = -1)
        @test_throws ArgumentError sampled(; nsamples = 0)
        @test_throws ArgumentError SampledSimulation(pendulum, 0.0; nsamples = 4)
        @test_throws ArgumentError SampledSimulation(
            pendulum,
            tstep;
            nsamples = 4,
            nsystem = 0,
        )
        approx = SampledSimulation(pendulum, tstep; nsamples = 4)
        for name in (:approx_mode, :overapproximation_map, :approx, :time_step, :nsystem)
            @test_throws ArgumentError build_sampled_abstraction(
                approx;
                grids...,
                NamedTuple{(name,)}((nothing,))...,
            )
        end
        # It is an approximation of Dionysos, and the system to use is its RK4 map.
        @test approx isa ST.DiscreteTimeSystemOverApproximation
        @test ST.get_system(approx) isa MS.ConstrainedBlackBoxControlDiscreteSystem
        image = ST.reach_set(
            approx,
            LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.3, 0.5)),
            SVector(0.0),
            nothing,
        )
        @test image isa SampledImage && length(image.points) == 4
    end

    @testset "double pendulum" begin
        double = DoublePendulum.system()
        dgrids = (;
            state_grid = MP.GridFree(
                SVector(0.0, 0.0, 0.0, 0.0),
                SVector(0.5, 0.5, 1.0, 1.0),
            ),
            input_grid = MP.GridFree(SVector(0.0), SVector(11.0)),
        )
        dtstep = 0.05
        dsampled(; kwargs...) =
            build_sampled_abstraction(double; dgrids..., time_step = dtstep, kwargs...)
        elapsed = @elapsed a = dsampled(; nsamples = 8, seed = 1)
        @test !isempty(transitions(a))
        @test elapsed < 120                                    # seconds, compilation included
        @test SY.get_n_state(abstract_system(a)) > 5000
        @test transitions(a) == transitions(dsampled(; nsamples = 8, seed = 1))
        @test transitions(a) == transitions(
            dsampled(; nsamples = 8, seed = 1, execution_backend = SY.ThreadedBackend(0.2)),
        )
        @test transitions(a) != transitions(dsampled(; nsamples = 8, seed = 2))
        center = BA.build_abstraction(
            double;
            dgrids...,
            time_step = dtstep,
            approx_mode = UGA.CENTER_SIMULATION,
        )
        @test transitions(dsampled(; nsamples = 0, include_center = true)) ==
              transitions(center)
        @test inside(
            successors(center),
            successors(dsampled(; nsamples = 8, include_center = true)),
        )
        @test inside(
            successors(dsampled(; nsamples = 2)),
            successors(dsampled(; nsamples = 16)),
        )

        # Four states: 17 points are the center and the 16 corners of every cell.
        corners(; kwargs...) = dsampled(; nsamples = 17, sampling = :uniform, kwargs...)
        @test transitions(dsampled(; nsamples = 1, sampling = :uniform)) ==
              transitions(center)
        @test inside(successors(center), successors(corners()))
        @test transitions(corners()) ==
              transitions(corners(; execution_backend = SY.ThreadedBackend(0.2)))
        @test inside(
            successors(dsampled(; nsamples = 5, sampling = :uniform)),
            successors(corners()),
        )
        # By hand, on a sample of the cells: the center and the corners, simulated with the RK4 map.
        opt = corners()
        sym = abstract_system(opt)
        found = successors(opt)
        F = ST.get_system_map(SampledSimulation(double, dtstep; nsamples = 1))
        mismatches = 0
        for s in 1:97:SY.get_n_state(sym), i in 1:SY.get_n_input(sym)
            cell, u = SY.get_concrete_elem(sym, s), SY.get_concrete_input(sym, i)
            c, ρ = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
            xs = vcat(
                [c],
                vec([
                    c .+ ρ .* SVector(s1, s2, s3, s4) for
                    s1 in (-1, 1), s2 in (-1, 1), s3 in (-1, 1), s4 in (-1, 1)
                ]),
            )
            @test length(xs) == 17
            targets = [SY.get_abstract_state(sym, F(x, u)) for x in xs]
            expected =
                all(q -> SY.is_allowed_state(sym, q), targets) ? Set(targets) : nothing
            mismatches += get(found, (s, i), nothing) == expected ? 0 : 1
        end
        @test mismatches == 0
    end
end
