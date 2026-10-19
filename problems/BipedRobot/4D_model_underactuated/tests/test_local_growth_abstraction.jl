# Tests for LocalGrowthAbstraction.jl. Run with an environment that has Dionysos, Symbolics,
# IntervalArithmetic, MathematicalSystems, LazySets, StaticArrays, JLD2 and MathOptInterface:
#   julia --project=<env> tests/test_local_growth_abstraction.jl
using Test
using Random
using StaticArrays
using Symbolics
using IntervalArithmetic          # loads Dionysos' interval/Symbolics extension
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
        opt_global =
            build_abstraction(lgb.global_growth; state_grid = grid, input_grid = igrid)
        sym_l = MOI.get(opt_local, MOI.RawOptimizerAttribute("abstract_system"))
        sym_g = MOI.get(opt_global, MOI.RawOptimizerAttribute("abstract_system"))

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
