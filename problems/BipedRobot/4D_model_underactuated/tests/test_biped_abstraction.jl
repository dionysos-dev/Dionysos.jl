# Tests for BipedAbstraction.jl (and CertifiedRK4.jl). The cache files and the chunks are also tried
# with an uncertified growth bound of the pendulum (`LocalGrowthBound.cellwise_bound`). Run with an
# environment that has Dionysos, Symbolics, IntervalArithmetic, CoreMath, MathematicalSystems, LazySets,
# StaticArrays, JLD2 and MathOptInterface:
#   julia --project=<env> tests/test_biped_abstraction.jl
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
import IntervalArithmetic

const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MP = Dionysos.Mapping
const PR = Dionysos.Problem
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction

include(joinpath(@__DIR__, "..", "LocalGrowthBound.jl"))
include(joinpath(@__DIR__, "..", "BipedAbstraction.jl"))
using .LocalGrowthBound
using .BipedAbstraction
const CertifiedRK4 = BipedAbstraction.CertifiedRK4

Random.seed!(1234)

# ---------------------------------------------------------------------------------------------
# A damped pendulum, ẋ = (ω, -a sin θ - b ω + u)
# ---------------------------------------------------------------------------------------------
@variables x[1:2] u[1:1]
f_pend = [x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u[1]]
X_pend = LazySets.Hyperrectangle(; low = [-3.0, -6.0], high = [3.0, 6.0])
U_pend = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
tstep = 0.1

dyn = SymbolicDynamics(f_pend, collect(x), collect(u))
system_pend = MS.ConstrainedBlackBoxControlContinuousSystem(dyn.f, 2, 1, X_pend, U_pend)

# The abstraction of the pendulum by the `GROWTH` mode of the Dionysos optimizer.
pendulum_abstraction(jacobian_bound; kwargs...) = build_abstraction(
    system_pend;
    approx_mode = UGA.GROWTH,
    jacobian_bound,
    time_step = tstep,
    kwargs...,
)

@testset "BipedAbstraction" begin
    @testset "certified RK4 successors" begin
        @test CertifiedRK4GrowthBound ===
              BipedAbstraction.CertifiedRK4.CertifiedRK4GrowthBound
        @variables z v
        X = LazySets.Hyperrectangle([0.0], [3.0])
        U = LazySets.Hyperrectangle([0.0], [1.0])
        ugrid = MP.GridFree(SVector(0.0), SVector(1.0))
        function check_successor(approx, grid, x, u; implicit = false)
            opt = build_abstraction(
                approx;
                state_grid = grid,
                input_grid = ugrid,
                use_implicit_mapping = implicit,
            )
            sym = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
            q = SY.get_abstract_state(sym, x)
            i = SY.get_abstract_input(sym, u)
            y = ST.get_system_map(approx)(x, u)
            target = SY.get_abstract_state(sym, y)
            successors = SY.post(sym, q, i)
            @test target !== nothing && target in successors
            return successors
        end

        # The ties-to-even quantizer maps 0.5 to index 0 and 1.5 to index 2.
        # The old deflated box for this action listed only index 1.
        translation = CertifiedRK4GrowthBound([v], [z], [v], X, U, 1.0; nsystem = 1)
        grid = MP.GridFree(SVector(0.0), SVector(1.0))
        @test length(check_successor(translation, grid, SVector(0.5), SVector(1.0))) == 3
        check_successor(translation, grid, SVector(0.5), SVector(1.0); implicit = true)
        # An odd source index excludes its tie faces. Its rounded updates may
        # still land on a tie, but they cannot reach the preceding odd cell.
        @test length(check_successor(translation, grid, SVector(1.0), SVector(1.0))) == 1

        # A contracting ODE can have an expanding RK4 step. This point is
        # inside the source cell and goes to a different index.
        nonlinear = CertifiedRK4GrowthBound([-tanh(10z)], [z], [v], X, U, 1.0; nsystem = 1)
        check_successor(
            nonlinear,
            MP.GridFree(SVector(0.0), SVector(0.2)),
            SVector(0.08),
            SVector(0.0),
        )

        # Powers. Symbolics hoists a literal exponent into a variable bound once
        # (`var"##cse#1" = 2`): it must still count as a constant exponent, and `z * z`
        # is rewritten to `z^2` by Symbolics itself. The Float64 map must agree with
        # an ordinary RK4 of the same field, and every successor must be covered.
        for (field, plain, x0, h) in (
            ([-z^3], (x, u) -> SVector(-x[1]^3), 0.8, 0.2),
            ([z^2 - 2 * z], (x, u) -> SVector(x[1]^2 - 2 * x[1]), 0.7, 0.1),
            ([(1 + z^2)^0.5], (x, u) -> SVector(sqrt(1 + x[1]^2)), 1.0, 0.25),
            ([-z * z], (x, u) -> SVector(-x[1]^2), 0.6, 0.05),
            (
                [1 / (z^2 + 0.5)^2 - z],
                (x, u) -> SVector(1 / (x[1]^2 + 0.5)^2 - x[1]),
                -0.4,
                0.1,
            ),
        )
            approx = CertifiedRK4GrowthBound(field, [z], [v], X, U, 0.5)
            @test only(ST.get_system_map(approx)(SVector(x0), SVector(0.0))) ≈ only(
                ST.runge_kutta4(
                    plain,
                    SVector(x0),
                    SVector(0.0),
                    0.5,
                    ST.DEFAULT_NUM_SUBSTEPS,
                ),
            ) rtol = 1e-12
            check_successor(
                approx,
                MP.GridFree(SVector(0.0), SVector(h)),
                SVector(x0),
                SVector(0.0),
            )
        end
        # A variable exponent has no sound interval image here.
        @test_throws ArgumentError CertifiedRK4GrowthBound([2^z], [z], [v], X, U, 0.1)

        # The local derivative should recover contraction that the direct
        # interval evaluation loses to dependency on the starting state.
        contraction = CertifiedRK4GrowthBound([-z], [z], [v], X, U, 0.5; nsystem = 1)
        img = ST.get_over_approximation_map(contraction)(
            LazySets.Hyperrectangle([0.0], [0.1]),
            SVector(0.0),
        )
        @test 0.06 < only(LazySets.high(img)) < 0.061
        @test -0.061 < only(LazySets.low(img)) < -0.06

        # The previous biped miss comes from a Float64 update of this form.
        translate_small = CertifiedRK4GrowthBound(
            [v],
            [z],
            [v],
            X,
            LazySets.Hyperrectangle([3.0], [0.0]),
            1 / 60,
        )
        grid = MP.GridFree(SVector(0.0), SVector(0.05))
        point = SVector(0.07500000000000001)
        cmd = SVector(3.0)
        opt = build_abstraction(
            translate_small;
            state_grid = grid,
            input_grid = MP.GridFree(SVector(3.0), SVector(3.0)),
        )
        sym = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
        y = ST.get_system_map(translate_small)(point, cmd)
        @test MP.get_pos_by_coord(grid, y) == (2,)
        @test SY.get_abstract_state(sym, y) in SY.post(
            sym,
            SY.get_abstract_state(sym, point),
            SY.get_abstract_input(sym, cmd),
        )

        # A partially defined stage cannot certify an action for the whole cell.
        partial = CertifiedRK4GrowthBound([sqrt(z)], [z], [v], X, U, 0.1)
        @test ST.get_over_approximation_map(partial)(
            LazySets.Hyperrectangle([0.0], [0.1]),
            SVector(0.0),
        ) === nothing
        opt = build_abstraction(partial; state_grid = grid, input_grid = ugrid)
        sym = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
        @test isempty(
            SY.post(
                sym,
                SY.get_abstract_state(sym, [0.0]),
                SY.get_abstract_input(sym, [0.0]),
            ),
        )
    end

    @testset "saving and loading the transitions" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        info = (; tstep, h = [0.4, 0.8])
        file = joinpath(mktempdir(), "abstraction.jld2")
        pairs_of(opt) = Set(
            SY.enum_transitions(MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))),
        )
        # Counts the enclosures computed, to see whether a build computes anything.
        cb = cellwise_bound(dyn)
        enclosures = Ref(0)
        counted = ST.CellwiseBound(
            (cell, uu, t) -> (enclosures[] += 1; cb.enclosure(cell, uu, t)),
            cb.bound,
        )

        computed = pendulum_abstraction(
            counted;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test isfile(file)                                     # written after the computation
        calls = enclosures[]
        @test calls > 0
        loaded = pendulum_abstraction(
            counted;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test enclosures[] == calls                            # nothing was computed
        @test pairs_of(loaded) == pairs_of(computed)
        @test !isempty(pairs_of(loaded))

        # The loaded abstraction is cached in the optimizer, ready for a control problem.
        @test MOI.get(loaded, MOI.RawOptimizerAttribute("abstract_system")) !== nothing

        # A file computed for other parameters, or another grid, is refused.
        @test_throws ErrorException pendulum_abstraction(
            counted;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = (; tstep, h = [0.4, 0.4]),
        )
        @test_throws ErrorException pendulum_abstraction(
            counted;
            state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.8)),
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
    end

    @testset "files of the first format (single arrays) still load" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        cb = cellwise_bound(dyn)
        full = pendulum_abstraction(cb; state_grid = grid, input_grid = igrid)
        sym = MOI.get(full, MOI.RawOptimizerAttribute("abstract_system"))
        tr = collect(SY.enum_transitions(sym))
        file = joinpath(mktempdir(), "legacy.jld2")
        JLD2.jldsave(
            file;
            n_state = SY.get_n_state(sym),
            n_input = SY.get_n_input(sym),
            info = (; tstep),
            probe = BipedAbstraction._id_probe(sym),
            target = Int32[t[1] for t in tr],
            source = Int32[t[2] for t in tr],
            symbol = Int32[t[3] for t in tr],
        )
        loaded = pendulum_abstraction(
            cb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = (; tstep),
        )
        sym2 = MOI.get(loaded, MOI.RawOptimizerAttribute("abstract_system"))
        @test Set(SY.enum_transitions(sym2)) == Set(tr)
    end

    @testset "chunks of a Slurm array, merged" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        cb = cellwise_bound(dyn)
        dir = mktempdir()
        nchunks = 3
        transitions(opt) = Set(
            SY.enum_transitions(MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))),
        )
        full = pendulum_abstraction(cb; state_grid = grid, input_grid = igrid)

        # One run per task of the array: each writes its chunk and caches nothing.
        for id in 1:nchunks
            pendulum_abstraction(
                cb;
                state_grid = grid,
                input_grid = igrid,
                execution_backend = SY.SlurmArrayBackend(
                    nchunks,
                    id,
                    dir,
                    :contiguous,
                    true,
                ),
            )
        end
        @test length(readdir(dir)) == nchunks

        # Without a file, the chunks are merged into the optimizer.
        in_memory = pendulum_abstraction(
            cb;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
        )
        @test transitions(in_memory) == transitions(full)
        @test !isempty(transitions(in_memory))

        # With a file, they go straight to it, without building the abstraction in memory.
        file = joinpath(dir, "merged.jld2")
        info = (; tstep)
        streamed = pendulum_abstraction(
            cb;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
            abstraction_file = file,
            abstraction_info = info,
        )
        @test isfile(file) && !isfile(file * ".tmp")
        @test MOI.get(streamed, MOI.RawOptimizerAttribute("abstract_system")) === nothing
        JLD2.jldopen(file, "r") do f
            @test f["nparts"] == nchunks
            @test f["n_transitions"] == length(transitions(full))
        end
        loaded = pendulum_abstraction(                            # ... which loads like a computation
            cb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test transitions(loaded) == transitions(full)

        # A chunk that is missing, or computed for another number of chunks, is an error.
        @test_throws ErrorException pendulum_abstraction(
            cb;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks = nchunks + 1),
        )
    end
    @testset "certified abstraction: cache files and Slurm chunks carry their provenance" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        info = (; tstep, h = [0.4, 0.8])
        cert =
            CertifiedRK4GrowthBound(f_pend, collect(x), collect(u), X_pend, U_pend, tstep)
        transitions(opt) = Set(
            SY.enum_transitions(MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))),
        )
        full = build_abstraction(cert; state_grid = grid, input_grid = igrid)
        @test !isempty(transitions(full))
        message(f) =
            try
                f()
                ""
            catch e
                sprint(showerror, e)
            end

        # What the provenance identifies.
        prov = BipedAbstraction.CertifiedRK4.get_provenance(cert)
        @test prov.approximation == "CertifiedRK4GrowthBound"
        @test prov.tstep == tstep && prov.nsystem == cert.nsystem
        @test prov.kernel_version == BipedAbstraction.CertifiedRK4.KERNEL_VERSION
        @test BipedAbstraction.provenance(cellwise_bound(dyn)) === nothing   # a growth bound of Dionysos keeps no record
        @test BipedAbstraction.CertifiedRK4.get_provenance(
            CertifiedRK4GrowthBound(f_pend, collect(x), collect(u), X_pend, U_pend, tstep),
        ) == prov                                          # the same dynamics, the same fingerprint
        f_other = [x[2], -9.81 * sin(x[1]) - 0.4 * x[2] + u[1]]
        other_dynamics =
            CertifiedRK4GrowthBound(f_other, collect(x), collect(u), X_pend, U_pend, tstep)
        other_tstep = CertifiedRK4GrowthBound(
            f_pend,
            collect(x),
            collect(u),
            X_pend,
            U_pend,
            tstep / 2,
        )
        @test BipedAbstraction.CertifiedRK4.get_provenance(other_dynamics).dynamics !=
              prov.dynamics
        @test BipedAbstraction.CertifiedRK4.get_provenance(other_tstep).tstep != prov.tstep
        @test BipedAbstraction.abstraction_provenance(
            MOI.get(full, MOI.RawOptimizerAttribute("abstract_system")),
            cellwise_bound(dyn),
        ) === nothing

        # A cache file: written with the provenance, loaded back, equal to the computation.
        file = joinpath(mktempdir(), "certified.jld2")
        computed = build_abstraction(
            cert;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test transitions(computed) == transitions(full)
        JLD2.jldopen(file, "r") do f
            @test f["provenance"].approximation == prov
        end
        loaded = build_abstraction(
            cert;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test transitions(loaded) == transitions(full)

        # An uncertified file, even with the same `info`, is refused (and the other way round).
        legacy_file = joinpath(mktempdir(), "uncertified.jld2")
        pendulum_abstraction(
            cellwise_bound(dyn);
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = legacy_file,
            abstraction_info = info,
        )
        JLD2.jldopen(legacy_file, "r") do f
            @test !haskey(f, "provenance")
        end
        refused = message() do
            build_abstraction(
                cert;
                state_grid = grid,
                input_grid = igrid,
                abstraction_file = legacy_file,
                abstraction_info = info,
            )
        end
        @test occursin("no provenance record", refused)
        refused = message() do
            pendulum_abstraction(
                cellwise_bound(dyn);
                state_grid = grid,
                input_grid = igrid,
                abstraction_file = file,
                abstraction_info = info,
            )
        end
        @test occursin("provenance of a certified construction", refused)

        # Other dynamics, another sampling period, another grid with the same `info`: refused.
        for (approx, g) in (
            (other_dynamics, grid),
            (other_tstep, grid),
            (cert, MP.GridFree(SVector(0.2, 0.0), SVector(0.4, 0.8))),
        )
            refused = message() do
                build_abstraction(
                    approx;
                    state_grid = g,
                    input_grid = igrid,
                    abstraction_file = file,
                    abstraction_info = info,
                )
            end
            @test occursin("another construction", refused)
        end

        # The chunks of a Slurm array: merged into the optimizer or into a file, equal to the computation.
        nchunks = 3
        # `approx` is a certified approximation, or a `CellwiseBound` (an uncertified growth bound).
        build_with(approx; kwargs...) =
            if approx isa CertifiedRK4GrowthBound
                build_abstraction(approx; state_grid = grid, input_grid = igrid, kwargs...)
            else
                pendulum_abstraction(
                    approx;
                    state_grid = grid,
                    input_grid = igrid,
                    kwargs...,
                )
            end
        function write_chunks(approx, dir)
            for id in 1:nchunks
                build_with(
                    approx;
                    execution_backend = SY.SlurmArrayBackend(
                        nchunks,
                        id,
                        dir,
                        :contiguous,
                        true,
                    ),
                )
            end
            return dir
        end
        dir = write_chunks(cert, mktempdir())
        in_memory = build_abstraction(
            cert;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
        )
        @test transitions(in_memory) == transitions(full)
        merged = joinpath(dir, "merged.jld2")
        streamed = build_abstraction(
            cert;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
            abstraction_file = merged,
            abstraction_info = info,
        )
        @test MOI.get(streamed, MOI.RawOptimizerAttribute("abstract_system")) === nothing
        JLD2.jldopen(merged, "r") do f
            @test f["provenance"].approximation == prov
        end
        @test transitions(
            build_abstraction(
                cert;
                state_grid = grid,
                input_grid = igrid,
                abstraction_file = merged,
                abstraction_info = info,
            ),
        ) == transitions(full)

        # Certified chunks are refused by an uncertified merge, and uncertified chunks by a
        # certified one; a chunk of other dynamics cannot slip into a set of certified chunks.
        @test_throws ErrorException pendulum_abstraction(
            cellwise_bound(dyn);
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
        )
        legacy_dir = write_chunks(cellwise_bound(dyn), mktempdir())
        refused = message() do
            build_abstraction(
                cert;
                state_grid = grid,
                input_grid = igrid,
                merge_chunks = (; outdir = legacy_dir, nchunks),
            )
        end
        @test occursin("no provenance record", refused)
        other_dir = write_chunks(other_dynamics, mktempdir())
        merge_from(outdir) = message() do
            build_abstraction(
                cert;
                state_grid = grid,
                input_grid = igrid,
                merge_chunks = (; outdir, nchunks),
            )
        end
        # A chunk file of other dynamics, copied over a certified chunk: it no longer is the file
        # that the record of the certified chunk describes.
        mixed = write_chunks(cert, mktempdir())
        cp(
            SY.slurm_chunk_file(other_dir, 2, nchunks),
            SY.slurm_chunk_file(mixed, 2, nchunks);
            force = true,
        )
        refused = merge_from(mixed)
        @test occursin("chunk 2 of $nchunks was modified", refused)
        # ... and with its own record, it is a chunk of another construction.
        cp(
            BipedAbstraction._provenance_file(other_dir, 2, nchunks),
            BipedAbstraction._provenance_file(mixed, 2, nchunks);
            force = true,
        )
        refused = merge_from(mixed)
        @test occursin("chunk 2 of $nchunks cannot be used", refused) &&
              occursin("another construction", refused)
        # Whoever recomputes a chunk without a record (here a cell-wise growth bound) removes the
        # record of the certified chunk it overwrites: the chunk cannot pass for a certified one.
        recomputed = write_chunks(cert, mktempdir())
        pendulum_abstraction(
            cellwise_bound(dyn);
            state_grid = grid,
            input_grid = igrid,
            execution_backend = SY.SlurmArrayBackend(
                nchunks,
                1,
                recomputed,
                :contiguous,
                true,
            ),
        )
        @test !isfile(BipedAbstraction._provenance_file(recomputed, 1, nchunks))
        @test occursin("no provenance record", merge_from(recomputed))
        # An untouched set of chunks still merges.
        @test merge_from(dir) == ""

        # Overriding what the certificate is about stays an error.
        @test_throws ArgumentError build_abstraction(
            cert;
            state_grid = grid,
            input_grid = igrid,
            approx_mode = nothing,
        )
    end

    @testset "fingerprints are stable, and provenance is checked when read back" begin
        fp = CertifiedRK4
        hex(x) = fp.fingerprint_string(fp.fingerprint(x))
        # Reference values from an independent implementation (Python): fingerprints are stored in
        # cache files, so these must never change (not with the Julia version either).
        @test hex("") == "985b2cc3d2245173"
        @test hex("abc") == "36a468722b4f6810"
        @test hex(42) == "b960a184f07032c6"
        @test hex(-1) == "685cd83ad34b3424"
        @test hex(1.5) == "0dcdddf54e95fb20"
        @test hex((1, 2, 3)) == "8c5e000ead7d66ca"
        @test hex(["x", 1, 2.5]) == "3528d36edbda2e4a"
        @test hex(("ab", "c")) == "710c416caa27ba5b"
        @test hex(("a", "bc")) == "3d3beb0ed0d7ace9"
        @test fp.fingerprint_string(fp.fingerprint(2, fp.fingerprint(1))) ==
              "a14c003dc26df266"
        @test fp.fingerprint(1) != fp.fingerprint(1.0) &&
              fp.fingerprint(1) != fp.fingerprint("1")
        @test fp.fingerprint(Int32(7)) == fp.fingerprint(7)
        @test fp.fingerprint(0.0) != fp.fingerprint(-0.0)
        file = joinpath(mktempdir(), "bytes")
        write(file, "abc")
        @test fp.fingerprint_file(file) == "e71fa2190541574b"      # plain FNV-1a 64 of the bytes "abc"

        recorded =
            (; approximation = (; approximation = "recorded", version = 1), ids = "abc")
        @test BipedAbstraction.check_provenance(recorded, recorded, "a file") === nothing
        @test BipedAbstraction.check_provenance(nothing, nothing, "a file") === nothing
        message(found, expected) =
            try
                BipedAbstraction.check_provenance(found, expected, "chunk 2")
                ""
            catch e
                sprint(showerror, e)
            end
        err = message(nothing, recorded)
        @test occursin("chunk 2 cannot be used", err) &&
              occursin("no provenance record", err)
        @test occursin("provenance of a certified construction", message(recorded, nothing))
        err = message((; recorded..., ids = "def"), recorded)
        @test occursin("another construction", err) && occursin("stored:", err)
    end
end
