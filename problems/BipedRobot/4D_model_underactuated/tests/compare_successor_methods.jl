# Comparison of four ways of computing the successors of a cell, on the simple and the double pendulum:
#
# * the global ("uniform") growth bound: one bound of the Jacobian for the whole state space (the one
#   of `problems/Pendulum`);
# * the cell-wise (local) growth bound, `LocalGrowthBound.cellwise_bound`;
# * `N` points of every cell chosen at random, simulated (`SampledAbstraction`, `sampling = :random`);
# * `N` points of every cell on a regular design (center, corners, ...), simulated (`sampling = :uniform`).
#
# The two sampled methods are run for three values of `N` per system (`study.Ns`: the middle one,
# `study.N`, is the reference), so that every system, sampling time and cell size gives eight
# abstractions. They are built with the Dionysos optimizer, and the script records the total and the
# average (per cell and input) number of transitions, the total and the average computation time,
# the share of (cell, input) pairs that keep a transition, and the histogram of the number of
# successors per pair. It also times the one-off symbolic setup of the cell-wise bound. The results go
# to `successor_methods_results/` (a Julia serialization and a CSV); `plot_successor_methods.jl` turns
# them into the figures and tables of the presentation.
#
# Timing: one thread, the build of the abstraction only (without the compilation: every method is run
# once on a small problem first, and a run shorter than a few seconds is repeated and the fastest taken).
# Run it with an environment that has Dionysos, Symbolics, IntervalArithmetic, CoreMath, JLD2, LazySets,
# MathematicalSystems, StaticArrays, MathOptInterface and Plots (the pendulum problems load it):
#
#   julia -t 1 --project=<env> tests/compare_successor_methods.jl

using Dates
using Printf
using Serialization
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
const BA = BipedAbstraction

const METHODS = (:global, :cellwise, :random, :uniform)
const SAMPLED = (:random, :uniform)

# The studies: the system, its constant bound of the Jacobian, the numbers `Ns` of points of the sampled
# methods (`N`: the reference one, the middle one), and the sampling times and cell sizes to try. The
# double pendulum is taken on the whole period of both angles, [-π, π]², and on |ω| ≤ 4 (the Jacobian
# bound of the problem file holds on |ω| ≤ 5, hence on this box). The angles are not wrapped: a successor
# beyond ±π is outside the domain and is dropped, as for every other state-space boundary.
#
# `Ns`: the uniform design is the center, then the 2ⁿ corners, then the 3ⁿ lattice, then the 5ⁿ one.
# Simple pendulum: 5 (center and corners), 9 (3 × 3 lattice, the reference), 25 (5 × 5 lattice).
# Double pendulum: 9 (center and half of the 16 corners), 17 (center and corners, the reference), 81 (the
# 3⁴ lattice; the 5⁴ lattice, 625 points, is out of reach for the largest grid).
function studies()
    simple = (
        name = :simple,
        system = SimplePendulum.system(),
        jacobian_bound = SimplePendulum.jacobian_bound(),
        N = 9,
        Ns = (5, 9, 25),
        taus = [0.05, 0.1, 0.2],
        sizes = [[0.2, 0.4], [0.1, 0.2], [0.05, 0.1]],
        input_step = 1.5,
        warm_size = [0.4, 0.8],
        warm_tau = 0.1,
    )
    box = LazySets.Hyperrectangle(;
        low = SVector(-π, -π, -4.0, -4.0),
        high = SVector(π, π, 4.0, 4.0),
    )
    double = (
        name = :double,
        system = DoublePendulum.system(; _X_ = box),
        jacobian_bound = DoublePendulum.jacobian_bound(),
        N = 17,
        Ns = (9, 17, 81),
        taus = [0.01, 0.02, 0.05],
        sizes = [[0.4, 0.4, 0.8, 0.8], [0.3, 0.3, 0.6, 0.6], [0.2, 0.2, 0.4, 0.4]],
        input_step = 11.0,
        warm_size = [0.8, 0.8, 1.6, 1.6],
        warm_tau = 0.01,
    )
    return [simple, double]
end

# The runs of a study, in the order of the figures: a variant is `(method, N)`, with `N = 0` for the
# growth bounds.
variants(study) = vcat(
    [(:global, 0), (:cellwise, 0)],
    [(:random, N) for N in study.Ns],
    [(:uniform, N) for N in study.Ns],
)

grids_of(study, size) = (;
    state_grid = MP.GridFree(
        zeros(SVector{length(size), Float64}),
        SVector{length(size)}(size),
    ),
    input_grid = MP.GridFree(SVector(0.0), SVector(study.input_step)),
)

# The abstraction of one method: the optimizer, with the abstraction cached. `N`: the number of points
# of the sampled methods (the reference one of the study by default).
function build(method::Symbol, study, dynamics, grids, tau, N = study.N)
    system = study.system
    if method === :global
        return BA.build_abstraction(
            system;
            grids...,
            time_step = tau,
            approx_mode = UGA.GROWTH,
            jacobian_bound = study.jacobian_bound,
        )
    elseif method === :cellwise
        return BA.build_abstraction(
            system;
            grids...,
            time_step = tau,
            approx_mode = UGA.GROWTH,
            jacobian_bound = cellwise_bound(dynamics),
        )
    elseif method === :random
        return build_sampled_abstraction(
            system;
            grids...,
            time_step = tau,
            nsamples = N,
            seed = 1,
        )
    elseif method === :uniform
        return build_sampled_abstraction(
            system;
            grids...,
            time_step = tau,
            nsamples = N,
            sampling = :uniform,
        )
    end
    return error("unknown method $method")
end

abstraction_of(opt) = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))

# The number of successors of every (source, input) pair, and its histogram (index k + 1: the number of
# pairs with k successors, k = 0 for the pairs without transition).
function successor_counts(opt)
    sym = abstraction_of(opt)
    counts = Dict{Tuple{Int, Int}, Int}()
    for (_, source, input) in SY.enum_transitions(sym)
        counts[(source, input)] = get(counts, (source, input), 0) + 1
    end
    pairs = SY.get_n_state(sym) * SY.get_n_input(sym)
    hist = zeros(Int, maximum(values(counts); init = 0) + 1)
    hist[1] = pairs - length(counts)
    for k in values(counts)
        hist[k + 1] += 1
    end
    return counts, hist
end

# The shortest of a few runs; a run that takes long is not repeated.
function timed(f; long = 2.0)
    GC.gc()
    t = @elapsed opt = f()
    best = t
    if t < long
        for _ in 1:(t < 0.2 ? 4 : 2)
            GC.gc()
            best = min(best, @elapsed(f()))
        end
    end
    return opt, best
end

function run_study(study)
    records = NamedTuple[]
    system = study.system
    # The one-off symbolic setup of the cell-wise bound for this system: the symbolic expressions of the
    # dynamics and of their derivatives, compiled, and used once (in a session where the symbolic
    # machinery is compiled already: see `main`).
    local dynamics
    setup = @elapsed begin
        dynamics = SymbolicDynamics(system)
        r0 = SVector{LazySets.dim(system.X), Float64}(study.warm_size) ./ 2
        cellwise_bound(dynamics).enclosure(
            LazySets.Hyperrectangle(zero(r0), r0),
            zeros(SVector{MS.inputdim(system), Float64}),
            study.warm_tau,
        )
    end
    # Compile every method on a small problem.
    for (method, N) in variants(study)
        build(
            method,
            study,
            dynamics,
            grids_of(study, study.warm_size),
            study.warm_tau,
            N == 0 ? study.N : N,
        )
    end
    for tau in study.taus, size in study.sizes
        grids = grids_of(study, size)
        for (method, N) in variants(study)
            opt, time = timed(() -> build(method, study, dynamics, grids, tau, N))
            sym = abstraction_of(opt)
            counts, hist = successor_counts(opt)
            n_cells, n_inputs = SY.get_n_state(sym), SY.get_n_input(sym)
            pairs = n_cells * n_inputs
            record = (;
                system = study.name,
                method,
                tau,
                size = copy(size),
                n_cells,
                n_inputs,
                pairs,
                N,
                transitions = SY.get_n_transitions(sym),
                enabled_pairs = length(counts),
                time,
                setup = method === :cellwise ? setup : 0.0,
                hist,
            )
            push!(records, record)
            @printf(
                "%-7s tau=%-5g h=%-26s %-9s %3d %9d transitions  %7.2f per pair  %5.1f%% enabled  %9.4f s  %8.2f µs/pair\n",
                study.name,
                tau,
                string(size),
                method,
                N,
                record.transitions,
                record.transitions / pairs,
                100 * record.enabled_pairs / pairs,
                time,
                1e6 * time / pairs
            )
            flush(stdout)
        end
    end
    return records, setup
end

function main()
    out = joinpath(@__DIR__, "successor_methods_results")
    mkpath(out)
    records = NamedTuple[]
    setups = Dict{Symbol, Float64}()
    # The first use of the symbolic machinery in a Julia session compiles it: timed on a throwaway
    # system, and reported apart from the setup of each system.
    cold = @elapsed begin
        dummy = SymbolicDynamics(
            SimplePendulum.system(; params = SimplePendulum.Params(; g = 5.0)),
        )
        cellwise_bound(dummy).enclosure(
            LazySets.Hyperrectangle(SVector(0.0, 0.0), SVector(0.1, 0.1)),
            SVector(0.0),
            0.1,
        )
    end
    @printf("first use of the symbolic machinery in this session: %.2f s\n", cold)
    for study in studies()
        r, setup = run_study(study)
        append!(records, r)
        setups[study.name] = setup
        @printf(
            "%s: one-off symbolic setup of the cell-wise bound: %.3f s\n",
            study.name,
            setup
        )
    end
    info = (;
        date = string(Dates.now()),
        julia = string(VERSION),
        cpu = Sys.cpu_info()[1].model,
        threads = Threads.nthreads(),
        sampled_points = Dict(s.name => collect(s.Ns) for s in studies()),
        reference_points = Dict(s.name => s.N for s in studies()),
        setups,
        cold_symbolic = cold,
    )
    serialize(joinpath(out, "results.jls"), (; records, info))
    open(joinpath(out, "results.csv"), "w") do io
        println(
            io,
            "system,method,tau,size,n_cells,n_inputs,pairs,N,transitions,enabled_pairs,time_s,setup_s",
        )
        for r in records
            println(
                io,
                join(
                    [
                        r.system,
                        r.method,
                        r.tau,
                        "\"" * join(r.size, " ") * "\"",
                        r.n_cells,
                        r.n_inputs,
                        r.pairs,
                        r.N,
                        r.transitions,
                        r.enabled_pairs,
                        r.time,
                        r.setup,
                    ],
                    ",",
                ),
            )
        end
    end
    println("written to $out")
    return nothing
end

# Run when it is the script, not when `plot_successor_methods.jl` includes it for the definitions.
(abspath(PROGRAM_FILE) == (@__FILE__)) && main()
