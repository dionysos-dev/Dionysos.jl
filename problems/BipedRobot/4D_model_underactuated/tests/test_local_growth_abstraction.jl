# Tests for LocalGrowthAbstraction.jl. Run with an environment that has Dionysos, Symbolics,
# IntervalArithmetic, CoreMath, MathematicalSystems, LazySets, StaticArrays, JLD2 and MathOptInterface:
#   julia --project=<env> tests/test_local_growth_abstraction.jl
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

include(joinpath(@__DIR__, "..", "Multibody2D.jl"))
include(joinpath(@__DIR__, "..", "LocalGrowthAbstraction.jl"))
using .Multibody2D
using .LocalGrowthAbstraction
const LGA = LocalGrowthAbstraction
const CertifiedRK4 = LGA.CertifiedRK4

Random.seed!(1234)

# ---------------------------------------------------------------------------------------------
# A damped pendulum, ẋ = (ω, -a sin θ - b ω + u): strongly non-linear only where |cos θ| is big
# ---------------------------------------------------------------------------------------------
@variables x[1:2] u[1:1]
f_pend = [x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u[1]]
X_pend = LazySets.Hyperrectangle(; low = [-3.0, -6.0], high = [3.0, 6.0])
U_pend = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
tstep = 0.1

lgb = LocalGrowthBound(f_pend, collect(x), collect(u), X_pend, U_pend, tstep)
dyn = lgb.dynamics
flow(x0, u0, τ) = ST.runge_kutta4(dyn.f, SVector{2}(x0), u0, τ, 400)   # accurate reference

function random_cell(h)
    c = SVector(rand() * 5 - 2.5, rand() * 10 - 5)
    return LazySets.Hyperrectangle(c, SVector(h...))
end
random_point(cell) =
    LazySets.center(cell) .+ (2 * rand(2) .- 1) .* LazySets.radius_hyperrectangle(cell)

@testset "LocalGrowthAbstraction" begin
    @testset "global Jacobian bound agrees with Dionysos' compute_jacobian_bound" begin
        sys = MS.ConstrainedBlackBoxControlContinuousSystem(
            (x, u) -> SVector(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u[1]),
            2,
            1,
            X_pend,
            U_pend,
        )
        dionysos_bound = ST.compute_jacobian_bound(sys)          # INPUT_BOUND
        for uval in (-2.0, 0.0, 1.3)
            L = LGA.jacobian_bound(dyn, [-3.0, -6.0], [3.0, 6.0], [uval])
            @test L ≈ dionysos_bound([uval])
        end
        # A narrower box gives a (weakly) smaller bound; a box where sin is monotone, much smaller.
        Lbig = LGA.jacobian_bound(dyn, [-3.0, -6.0], [3.0, 6.0], [0.0])
        Lsmall = LGA.jacobian_bound(dyn, [1.4, -1.0], [1.7, 1.0], [0.0])
        @test all(Lsmall .<= Lbig)
        @test Lsmall[2, 1] < 0.2 * Lbig[2, 1]                      # |cos θ| ≈ 0 near π/2
        # The "model" of Dionysos can also be traced from the black box:
        @test LGA.jacobian_bound(SymbolicDynamics(sys), [-3.0, -6.0], [3.0, 6.0], [0.0]) ≈
              Lbig
    end

    @testset "flow enclosure contains the whole tube" begin
        for _ in 1:20
            cell = random_cell((0.2, 0.4))
            uval = [4 * rand() - 2]
            global_img =
                ST.get_over_approximation_map(lgb.global_growth)(cell, SVector(uval...))
            B = LGA.flow_enclosure(
                dyn,
                LazySets.low(cell),
                LazySets.high(cell),
                LazySets.low(global_img),
                LazySets.high(global_img),
                uval,
                tstep,
            )
            @test B !== nothing
            for _ in 1:10, t in range(0, tstep; length = 11)
                y = flow(random_point(cell), SVector(uval...), t)
                @test all(B[1] .- 1e-9 .<= y .<= B[2] .+ 1e-9)
            end
        end
    end

    @testset "flow enclosure validates rounding and the domain" begin
        @variables z v
        translation = SymbolicDynamics([v], [z], [v])
        # The exact upper endpoint is 1 + 2^-53. Round-to-nearest evaluates it to 1.
        @test flow_enclosure(
            translation,
            [0.5],
            [1.0],
            [0.0],
            [1.0],
            [2.0^-53],
            1.0;
            max_iterations = 1,
        ) === nothing
        lo, hi = flow_enclosure(
            translation,
            [0.5],
            [1.0],
            [0.0],
            [1.0],
            [2.0^-53],
            1.0;
            shrink = true,
        )
        @test lo == [0.5]
        @test hi == [nextfloat(1.0)]

        # The exact upper face of the stored center-radius cell is also 1 + 2^-53.
        cell = LazySets.Hyperrectangle([1.0], [2.0^-53])
        _, cell_hi =
            flow_enclosure(translation, cell, [0.0], [2.0], [0.0], 1.0; shrink = true)
        @test cell_hi == [nextfloat(1.0)]

        # sqrt([-1, 2]) is bounded on its defined part. The box cannot certify an ODE
        # with sqrt(z) from negative initial states.
        partial = SymbolicDynamics([sqrt(z)], [z], [v])
        @test flow_enclosure(partial, [-0.5], [0.5], [-1.0], [2.0], [0.0], 0.1) === nothing
    end

    @testset "certified RK4 successors" begin
        @test CertifiedRK4GrowthBound === LGA.CertifiedRK4.CertifiedRK4GrowthBound
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

    @testset "local image is sound, and never worse than the global one" begin
        reset_stats!(lgb)
        over_local = ST.get_over_approximation_map(lgb)
        over_global = ST.get_over_approximation_map(lgb.global_growth)
        n_better = 0
        for _ in 1:60
            cell = random_cell((0.2, 0.4))
            uval = SVector(4 * rand() - 2)
            Rl, Rg = over_local(cell, uval), over_global(cell, uval)
            @test all(
                LazySets.radius_hyperrectangle(Rl) .<=
                LazySets.radius_hyperrectangle(Rg) .+ 1e-12,
            )
            @test LazySets.center(Rl) ≈ LazySets.center(Rg)
            n_better += any(
                LazySets.radius_hyperrectangle(Rl) .<
                LazySets.radius_hyperrectangle(Rg) .- 1e-9,
            )
            for _ in 1:20                                      # soundness: sampled successors
                y = flow(random_point(cell), uval, tstep)
                @test all(LazySets.low(Rl) .- 1e-9 .<= y .<= LazySets.high(Rl) .+ 1e-9)
            end
        end
        @test n_better > 30
        s = stats(lgb)
        @test s.calls == 60 && s.improved == n_better && s.enclosure_failures == 0
    end

    @testset "refinements only tighten" begin
        lgb2 = LocalGrowthBound(
            f_pend,
            collect(x),
            collect(u),
            X_pend,
            U_pend,
            tstep;
            refinements = 3,
        )
        lgb0 = LocalGrowthBound(
            f_pend,
            collect(x),
            collect(u),
            X_pend,
            U_pend,
            tstep;
            refinements = 0,
        )
        for _ in 1:20
            cell, uval = random_cell((0.2, 0.4)), SVector(4 * rand() - 2)
            r0 = LazySets.radius_hyperrectangle(
                ST.get_over_approximation_map(lgb0)(cell, uval),
            )
            r1 = LazySets.radius_hyperrectangle(
                ST.get_over_approximation_map(lgb)(cell, uval),
            )
            r3 = LazySets.radius_hyperrectangle(
                ST.get_over_approximation_map(lgb2)(cell, uval),
            )
            @test all(r3 .<= r1 .+ 1e-12) && all(r1 .<= r0 .+ 1e-12)
            r_glob = LazySets.radius_hyperrectangle(
                ST.get_over_approximation_map(lgb.global_growth)(cell, uval),
            )
            @test r0 ≈ r_glob                                   # 0 refinements = Dionysos' global bound
        end
    end

    @testset "fallback to the global bound when no enclosure exists" begin
        slow = LocalGrowthBound(
            f_pend,
            collect(x),
            collect(u),
            X_pend,
            U_pend,
            3.0;
            max_enclosure_iterations = 2,
        )
        cell, uval = random_cell((0.2, 0.4)), SVector(1.0)
        Rl = ST.get_over_approximation_map(slow)(cell, uval)
        Rg = ST.get_over_approximation_map(slow.global_growth)(cell, uval)
        @test LazySets.radius_hyperrectangle(Rl) ≈ LazySets.radius_hyperrectangle(Rg)
        @test stats(slow).enclosure_failures == 1
    end

    @testset "images touching a neighbour along a face do not overlap it" begin
        # ẋ = u on a lattice: τ·u is exactly one cell, so the image of a cell is exactly a cell.
        @variables y[1:2] w[1:2]
        Xl = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [1.0, 1.0])
        Ul = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [1.0, 1.0])
        h = SVector(0.1, 0.1)
        cell = LazySets.Hyperrectangle(SVector(0.0, 0.0), h ./ 2)
        grid = MP.GridFree(SVector(0.0, 0.0), h)
        ugrid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
        function nsuccessors(tol)
            l = LocalGrowthBound(
                collect(w),
                collect(y),
                collect(w),
                Xl,
                Ul,
                0.1;
                boundary_tolerance = tol,
            )
            sym = MOI.get(
                build_abstraction(l; state_grid = grid, input_grid = ugrid),
                MOI.RawOptimizerAttribute("abstract_system"),
            )
            q = SY.get_abstract_state(sym, [0.0, 0.0])
            i = SY.get_abstract_input(sym, [1.0, 0.0])
            return count(t -> t[2] == q && t[3] == i, SY.enum_transitions(sym))
        end
        @test nsuccessors(1e-9) == 1            # exactly the cell one step to the right
        @test nsuccessors(0.0) == 9             # closed cells: touching neighbours are included
    end

    @testset "transition systems: local ⊆ global, and sound" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        reset_stats!(lgb)
        opt_local = build_abstraction(lgb; state_grid = grid, input_grid = igrid)
        global_radius_calls = Ref(0)
        global_growth = ST.DiscreteTimeGrowthBound(
            ST.get_system(lgb.global_growth),
            (r, u) -> begin
                global_radius_calls[] += 1
                lgb.global_growth.growthbound_map(r, u)
            end,
        )
        opt_global = build_abstraction(global_growth; state_grid = grid, input_grid = igrid)
        sym_l = MOI.get(opt_local, MOI.RawOptimizerAttribute("abstract_system"))
        sym_g = MOI.get(opt_global, MOI.RawOptimizerAttribute("abstract_system"))
        # The optimizer wraps the function of `overapproximation_map` around the same system.
        @test ST.get_system(
            MOI.get(
                opt_local,
                MOI.RawOptimizerAttribute("discrete_time_system_approximation"),
            ),
        ) === ST.get_system(lgb)
        @test global_radius_calls[] == length(collect(SY.enum_inputs(sym_g)))

        targets(sym) = begin
            d = Dict{Tuple{Int, Int}, Set{Int}}()
            for (t, s, i) in SY.enum_transitions(sym)
                push!(get!(d, (s, i), Set{Int}()), t)
            end
            d
        end
        Tl, Tg = targets(sym_l), targets(sym_g)
        @test length(Tl) >= length(Tg)                         # every enabled pair stays enabled
        @test all(k -> haskey(Tl, k) && Tl[k] ⊆ Tg[k], keys(Tg))
        nl, ng = sum(length, values(Tl)), sum(length, values(Tg))
        avg_l, avg_g = nl / length(Tl), ng / length(Tg)
        @info "pendulum abstraction" cells = SY.get_n_state(sym_l) inputs =
            SY.get_n_input(sym_l) transitions_global = ng transitions_local = nl avg_successors_global =
            avg_g avg_successors_local = avg_l stats = stats(lgb)
        @test avg_l < avg_g                                    # fewer spurious transitions

        # Soundness: a sampled successor always lies in a target cell.
        for ((q, i), tg) in Tl
            x0 = SY.get_concrete_state(sym_l, q)
            uc = SY.get_concrete_input(sym_l, i)
            h = MP.get_h(MP.get_grid(SY.get_state_mapping(sym_l)))
            for _ in 1:3
                p = x0 .+ (2 * rand(2) .- 1) .* h ./ 2
                y = flow(p, uc, tstep)
                t = SY.get_abstract_state(sym_l, y)
                t === nothing || @test t in tg
            end
        end
    end

    @testset "saving and loading the transitions" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        info = (; tstep, h = [0.4, 0.8])
        file = joinpath(mktempdir(), "abstraction.jld2")
        pairs_of(opt) = Set(
            SY.enum_transitions(MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))),
        )

        computed = build_abstraction(
            lgb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test isfile(file)                                     # written after the computation
        calls = lgb.stats.calls[]
        loaded = build_abstraction(
            lgb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test lgb.stats.calls[] == calls                       # nothing was computed
        @test pairs_of(loaded) == pairs_of(computed)
        @test !isempty(pairs_of(loaded))

        # The loaded abstraction is cached in the optimizer, ready for a control problem.
        @test MOI.get(loaded, MOI.RawOptimizerAttribute("abstract_system")) !== nothing

        # A file computed for other parameters, or another grid, is refused.
        @test_throws ErrorException build_abstraction(
            lgb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = (; tstep, h = [0.4, 0.4]),
        )
        @test_throws ErrorException build_abstraction(
            lgb;
            state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.8)),
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
    end

    @testset "files of the first format (single arrays) still load" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        full = build_abstraction(lgb; state_grid = grid, input_grid = igrid)
        sym = MOI.get(full, MOI.RawOptimizerAttribute("abstract_system"))
        tr = collect(SY.enum_transitions(sym))
        file = joinpath(mktempdir(), "legacy.jld2")
        JLD2.jldsave(
            file;
            n_state = SY.get_n_state(sym),
            n_input = SY.get_n_input(sym),
            info = (; tstep),
            probe = LGA._id_probe(sym),
            target = Int32[t[1] for t in tr],
            source = Int32[t[2] for t in tr],
            symbol = Int32[t[3] for t in tr],
        )
        loaded = build_abstraction(
            lgb;
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
        dir = mktempdir()
        nchunks = 3
        transitions(opt) = Set(
            SY.enum_transitions(MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))),
        )
        full = build_abstraction(lgb; state_grid = grid, input_grid = igrid)

        # One run per task of the array: each writes its chunk and caches nothing.
        for id in 1:nchunks
            build_abstraction(
                lgb;
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
        in_memory = build_abstraction(
            lgb;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
        )
        @test transitions(in_memory) == transitions(full)
        @test !isempty(transitions(in_memory))

        # With a file, they go straight to it, without building the abstraction in memory.
        file = joinpath(dir, "merged.jld2")
        info = (; tstep)
        streamed = build_abstraction(
            lgb;
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
        loaded = build_abstraction(                            # ... which loads like a computation
            lgb;
            state_grid = grid,
            input_grid = igrid,
            abstraction_file = file,
            abstraction_info = info,
        )
        @test transitions(loaded) == transitions(full)

        # A chunk that is missing, or computed for another number of chunks, is an error.
        @test_throws ErrorException build_abstraction(
            lgb;
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
        prov = LGA.CertifiedRK4.get_provenance(cert)
        @test prov.approximation == "CertifiedRK4GrowthBound"
        @test prov.tstep == tstep && prov.nsystem == cert.nsystem
        @test prov.kernel_version == LGA.CertifiedRK4.KERNEL_VERSION
        @test LGA.provenance(lgb) === nothing           # the local growth bound keeps no record
        @test LGA.CertifiedRK4.get_provenance(
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
        @test LGA.CertifiedRK4.get_provenance(other_dynamics).dynamics != prov.dynamics
        @test LGA.CertifiedRK4.get_provenance(other_tstep).tstep != prov.tstep
        @test LGA.abstraction_provenance(
            MOI.get(full, MOI.RawOptimizerAttribute("abstract_system")),
            lgb,
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
        build_abstraction(
            lgb;
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
            build_abstraction(
                lgb;
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
        function write_chunks(approx, dir)
            for id in 1:nchunks
                build_abstraction(
                    approx;
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
        @test_throws ErrorException build_abstraction(
            lgb;
            state_grid = grid,
            input_grid = igrid,
            merge_chunks = (; outdir = dir, nchunks),
        )
        legacy_dir = write_chunks(lgb, mktempdir())
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
            LGA._provenance_file(other_dir, 2, nchunks),
            LGA._provenance_file(mixed, 2, nchunks);
            force = true,
        )
        refused = merge_from(mixed)
        @test occursin("chunk 2 of $nchunks cannot be used", refused) &&
              occursin("another construction", refused)
        # Whoever recomputes a chunk without a record (here the local growth bound) removes the
        # record of the certified chunk it overwrites: the chunk cannot pass for a certified one.
        recomputed = write_chunks(cert, mktempdir())
        build_abstraction(
            lgb;
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
        @test !isfile(LGA._provenance_file(recomputed, 1, nchunks))
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
        @test LGA.check_provenance(recorded, recorded, "a file") === nothing
        @test LGA.check_provenance(nothing, nothing, "a file") === nothing
        message(found, expected) =
            try
                LGA.check_provenance(found, expected, "chunk 2")
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

    @testset "acrobot normal form from Multibody2D" begin
        @variables q[1:2] qd[1:2] l1 l2 lc1 lc2 m1 m2 I1 I2 g
        sys = MultibodySystem(collect(q), collect(qd); gravity = g)
        link1 = add_frame!(sys, :link1, sys.world, [0, 0], q[1])
        link2 = add_frame!(sys, :link2, link1, [l1, 0], q[2])
        add_point!(sys, :body1, link1, [lc1, 0]; mass = m1, inertia = I1)
        add_point!(sys, :body2, link2, [lc2, 0]; mass = m2, inertia = I2)
        nf = normal_form(sys, selection_matrix(2, [2]))            # x = (q1, q2, q̇2, p₁), u = τ

        # Numerical values of the parameters are only substituted now.
        params = Dict(
            l1 => 1.0,
            l2 => 1.0,
            lc1 => 0.5,
            lc2 => 0.5,
            m1 => 1.0,
            m2 => 1.0,
            I1 => 0.083,
            I2 => 0.083,
            g => 9.81,
        )
        xdot = Symbolics.substitute.(nf.xdot, Ref(params))
        # Sampling period and cells are chosen so that the Picard enclosure exists: the normal
        # form has large gains (∂v̇₂/∂q₂ ≈ 600 over the box), so τ·‖L‖ must stay small.
        xr = [0.4, 0.4, 0.8, 0.8]
        X = LazySets.Hyperrectangle(; low = -xr, high = xr)
        U = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
        tstep_a = 0.01
        h = SVector(0.1, 0.1, 0.2, 0.2)
        lgb_a = LocalGrowthBound(xdot, nf.state, nf.input, X, U, tstep_a)

        # The local bound does not degrade when the state set grows, the global one does.
        X_big = LazySets.Hyperrectangle(; low = 2.5 .* -xr, high = 2.5 .* xr)
        lgb_big = LocalGrowthBound(xdot, nf.state, nf.input, X_big, U, tstep_a)
        cell0 = LazySets.Hyperrectangle(SVector(0.1, 0.1, 0.2, 0.2), h ./ 2)
        rad(l, c, uu) =
            LazySets.radius_hyperrectangle(ST.get_over_approximation_map(l)(c, uu))
        uu = SVector(2.0)
        @test rad(lgb_big.global_growth, cell0, uu)[3] >
              1.5 * rad(lgb_a.global_growth, cell0, uu)[3]
        @test rad(lgb_big, cell0, uu)[3] < 1.1 * rad(lgb_a, cell0, uu)[3]
        @test rad(lgb_a, cell0, uu)[3] < 0.7 * rad(lgb_a.global_growth, cell0, uu)[3]
        reset_stats!(lgb_a)

        grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), h)
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        sym_l = MOI.get(
            build_abstraction(lgb_a; state_grid = grid, input_grid = igrid),
            MOI.RawOptimizerAttribute("abstract_system"),
        )
        sym_g = MOI.get(
            build_abstraction(lgb_a.global_growth; state_grid = grid, input_grid = igrid),
            MOI.RawOptimizerAttribute("abstract_system"),
        )
        ntl, ntg = SY.get_n_transitions(sym_l), SY.get_n_transitions(sym_g)
        @info "acrobot abstraction" cells = SY.get_n_state(sym_l) transitions_global = ntg transitions_local =
            ntl stats = stats(lgb_a)
        @test ntg > 0 && ntl > 0
        @test stats(lgb_a).enclosure_failures == 0

        pairs(sym) = begin
            d = Dict{Tuple{Int, Int}, Set{Int}}()
            for (t, s, i) in SY.enum_transitions(sym)
                push!(get!(d, (s, i), Set{Int}()), t)
            end
            d
        end
        Pl, Pg = pairs(sym_l), pairs(sym_g)
        @test all(k -> haskey(Pl, k) && Pl[k] ⊆ Pg[k], keys(Pg))
        @test sum(length, values(Pl)) / length(Pl) < sum(length, values(Pg)) / length(Pg)

        # Soundness against an accurate simulation of the same normal form.
        fa = lgb_a.dynamics.f
        for ((q0, i), tg) in Iterators.take(shuffle(collect(Pl)), 300)
            x0 = SY.get_concrete_state(sym_l, q0)
            uc = SY.get_concrete_input(sym_l, i)
            p = x0 .+ (2 * rand(4) .- 1) .* h ./ 2
            y = ST.runge_kutta4(fa, SVector{4}(p), uc, tstep_a, 200)
            t = SY.get_abstract_state(sym_l, y)
            t === nothing || @test t in tg
        end
    end
end
