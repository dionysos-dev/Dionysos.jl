# Tests for LocalGrowthBound.jl. Run with an environment that has Dionysos, Symbolics, IntervalArithmetic,
# MathematicalSystems, LazySets, StaticArrays and MathOptInterface:
#   julia --project=<env> tests/test_local_growth_bound.jl
using Test
using Random
using StaticArrays
using Symbolics
using IntervalArithmetic          # loads Dionysos' interval/Symbolics extension
import Dionysos
import LazySets
import MathematicalSystems as MS
import MathOptInterface as MOI

const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MP = Dionysos.Mapping
const PR = Dionysos.Problem
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction

include(joinpath(@__DIR__, "..", "Multibody2D.jl"))
include(joinpath(@__DIR__, "..", "LocalGrowthBound.jl"))
using .Multibody2D
using .LocalGrowthBound

Random.seed!(1234)

# ---------------------------------------------------------------------------------------------
# A damped pendulum, ẋ = (ω, -a sin θ - b ω + u): strongly non-linear only where |cos θ| is big
# ---------------------------------------------------------------------------------------------
@variables x[1:2] u[1:1]
f_pend = [x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u[1]]
X_pend = LazySets.Hyperrectangle(; low = [-3.0, -6.0], high = [3.0, 6.0])
U_pend = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
tstep = 0.1

dyn = SymbolicDynamics(f_pend, collect(x), collect(u))
system_pend = MS.ConstrainedBlackBoxControlContinuousSystem(dyn.f, 2, 1, X_pend, U_pend)
flow(x0, u0, τ) = ST.runge_kutta4(dyn.f, SVector{2}(x0), u0, τ, 400)   # accurate reference

# The growth bound of Dionysos over the whole state set (one matrix per input) is the baseline of
# the growth bound per cell, `cellwise_bound`, that Dionysos integrates in the same way.
global_bound(dynamics, X) =
    u -> jacobian_bound_over_box(dynamics, LazySets.low(X), LazySets.high(X), u)
growth(jacobian_bound; system = system_pend, tstep = tstep) =
    ST.discretize(ST.ContinuousTimeGrowthBound(system; jacobian_bound), tstep)
global_growth = growth(global_bound(dyn, X_pend))
cellwise_growth = growth(cellwise_bound(dyn))

# The transition system of `system` by the optimizer of Dionysos (`GROWTH` mode, the growth bound
# is the `jacobian_bound` keyword).
function growth_abstraction(system; state_grid, input_grid, attributes...)
    optimizer = MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(system, system.X),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    for (name, value) in attributes
        MOI.set(optimizer, MOI.RawOptimizerAttribute(String(name)), value)
    end
    MOI.optimize!(optimizer)
    return optimizer
end

# The abstraction of the pendulum by the `GROWTH` mode of the Dionysos optimizer.
pendulum_abstraction(jacobian_bound; kwargs...) = growth_abstraction(
    system_pend;
    approx_mode = UGA.GROWTH,
    jacobian_bound,
    time_step = tstep,
    kwargs...,
)

function random_cell(h)
    c = SVector(rand() * 5 - 2.5, rand() * 10 - 5)
    return LazySets.Hyperrectangle(c, SVector(h...))
end
random_point(cell) =
    LazySets.center(cell) .+ (2 * rand(2) .- 1) .* LazySets.radius_hyperrectangle(cell)

@testset "LocalGrowthBound" begin
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
            L = jacobian_bound_over_box(dyn, [-3.0, -6.0], [3.0, 6.0], [uval])
            @test L ≈ dionysos_bound([uval])
        end
        # A narrower box gives a (weakly) smaller bound; a box where sin is monotone, much smaller.
        Lbig = jacobian_bound_over_box(dyn, [-3.0, -6.0], [3.0, 6.0], [0.0])
        Lsmall = jacobian_bound_over_box(dyn, [1.4, -1.0], [1.7, 1.0], [0.0])
        @test all(Lsmall .<= Lbig)
        @test Lsmall[2, 1] < 0.2 * Lbig[2, 1]                      # |cos θ| ≈ 0 near π/2
        # The "model" of Dionysos can also be traced from the black box:
        @test jacobian_bound_over_box(
            SymbolicDynamics(sys),
            [-3.0, -6.0],
            [3.0, 6.0],
            [0.0],
        ) ≈ Lbig
    end

    @testset "Picard enclosure contains the whole tube" begin
        for _ in 1:20
            cell = random_cell((0.2, 0.4))
            uval = [4 * rand() - 2]
            global_img =
                ST.get_over_approximation_map(global_growth)(cell, SVector(uval...))
            B = picard_enclosure(
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

    @testset "Picard enclosure validates rounding and the domain" begin
        @variables z v
        translation = SymbolicDynamics([v], [z], [v])
        # The exact upper endpoint is 1 + 2^-53. Round-to-nearest evaluates it to 1.
        @test picard_enclosure(
            translation,
            [0.5],
            [1.0],
            [0.0],
            [1.0],
            [2.0^-53],
            1.0;
            max_iterations = 1,
        ) === nothing
        lo, hi = picard_enclosure(
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
            picard_enclosure(translation, cell, [0.0], [2.0], [0.0], 1.0; shrink = true)
        @test cell_hi == [nextfloat(1.0)]

        # sqrt([-1, 2]) is bounded on its defined part. The box cannot certify an ODE
        # with sqrt(z) from negative initial states.
        partial = SymbolicDynamics([sqrt(z)], [z], [v])
        @test picard_enclosure(partial, [-0.5], [0.5], [-1.0], [2.0], [0.0], 0.1) ===
              nothing
    end

    @testset "cellwise_bound supplies the two functions of ST.CellwiseBound" begin
        cb = cellwise_bound(dyn)
        @test cb isa ST.CellwiseBound
        cell, uval = random_cell((0.2, 0.4)), SVector(0.7)
        E = cb.enclosure(cell, uval, tstep)
        @test E isa FlowEnclosure
        for _ in 1:10, t in range(0, tstep; length = 11)         # the whole tube, not only its end
            y = flow(random_point(cell), uval, t)
            @test all(E.low .- 1e-9 .<= y .<= E.high .+ 1e-9)
        end
        @test cb.bound(E, uval) == E.L == jacobian_bound_over_box(dyn, E.low, E.high, uval)

        # No verified enclosure (a step too long for the Picard test): `nothing`, which Dionysos
        # turns into an input that is not admissible at that cell.
        slow = cellwise_bound(dyn; max_iterations = 2)
        @test slow.enclosure(cell, SVector(1.0), 3.0) === nothing
        over_slow = ST.get_over_approximation_map(growth(slow; tstep = 3.0))
        @test over_slow(cell, SVector(1.0)) === nothing

        # An enclosure whose Jacobian is not bounded disables the input too.
        f_unbounded = (xv, uv) -> SVector(uv[1])
        f_interval = (xb, ub) -> [ub[1]]
        jacobian_interval =
            (xb, ub) -> reshape([IntervalArithmetic.entireinterval(Float64)], 1, 1)
        unbounded = SymbolicDynamics{
            1,
            1,
            typeof(f_unbounded),
            typeof(f_interval),
            typeof(jacobian_interval),
        }(
            f_unbounded,
            f_interval,
            jacobian_interval,
        )
        @test cellwise_bound(unbounded).enclosure(
            LazySets.Hyperrectangle([0.0], [0.1]),
            SVector(1.0),
            0.1,
        ) === nothing
    end

    @testset "cellwise image is sound, and never worse than the global one" begin
        over_local = ST.get_over_approximation_map(cellwise_growth)
        over_global = ST.get_over_approximation_map(global_growth)
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
    end

    @testset "transition systems: cellwise ⊆ global, and sound" begin
        grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.8))
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        opt_local =
            pendulum_abstraction(cellwise_bound(dyn); state_grid = grid, input_grid = igrid)
        opt_global = pendulum_abstraction(
            global_bound(dyn, X_pend);
            state_grid = grid,
            input_grid = igrid,
        )
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
            SY.get_n_input(sym_l) transitions_global = ng transitions_cellwise = nl avg_successors_global =
            avg_g avg_successors_cellwise = avg_l
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
        dyn_a = SymbolicDynamics(xdot, nf.state, nf.input)
        system_a(X) = MS.ConstrainedBlackBoxControlContinuousSystem(dyn_a.f, 4, 1, X, U)
        growth_a(jacobian_bound, X) =
            growth(jacobian_bound; system = system_a(X), tstep = tstep_a)

        # The cellwise bound does not degrade when the state set grows, the global one does.
        X_big = LazySets.Hyperrectangle(; low = 2.5 .* -xr, high = 2.5 .* xr)
        cell0 = LazySets.Hyperrectangle(SVector(0.1, 0.1, 0.2, 0.2), h ./ 2)
        rad(approx, c, uu) =
            LazySets.radius_hyperrectangle(ST.get_over_approximation_map(approx)(c, uu))
        uu = SVector(2.0)
        cellwise_a = growth_a(cellwise_bound(dyn_a), X)
        cellwise_big = growth_a(cellwise_bound(dyn_a), X_big)
        global_a = growth_a(global_bound(dyn_a, X), X)
        global_big = growth_a(global_bound(dyn_a, X_big), X_big)
        @test rad(global_big, cell0, uu)[3] > 1.5 * rad(global_a, cell0, uu)[3]
        @test rad(cellwise_big, cell0, uu)[3] < 1.1 * rad(cellwise_a, cell0, uu)[3]
        @test rad(cellwise_a, cell0, uu)[3] < 0.7 * rad(global_a, cell0, uu)[3]

        grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), h)
        igrid = MP.GridFree(SVector(0.0), SVector(2.0))
        abstraction(jacobian_bound) = MOI.get(
            growth_abstraction(
                system_a(X);
                approx_mode = UGA.GROWTH,
                jacobian_bound,
                time_step = tstep_a,
                state_grid = grid,
                input_grid = igrid,
            ),
            MOI.RawOptimizerAttribute("abstract_system"),
        )
        sym_l = abstraction(cellwise_bound(dyn_a))
        sym_g = abstraction(global_bound(dyn_a, X))
        ntl, ntg = SY.get_n_transitions(sym_l), SY.get_n_transitions(sym_g)
        @info "acrobot abstraction" cells = SY.get_n_state(sym_l) transitions_global = ntg transitions_cellwise =
            ntl
        @test ntg > 0 && ntl > 0

        pairs(sym) = begin
            d = Dict{Tuple{Int, Int}, Set{Int}}()
            for (t, s, i) in SY.enum_transitions(sym)
                push!(get!(d, (s, i), Set{Int}()), t)
            end
            d
        end
        Pl, Pg = pairs(sym_l), pairs(sym_g)
        # No pair is disabled for lack of an enclosure, and every successor set is inside the global one.
        @test all(k -> haskey(Pl, k) && Pl[k] ⊆ Pg[k], keys(Pg))
        @test sum(length, values(Pl)) / length(Pl) < sum(length, values(Pg)) / length(Pg)

        # Soundness against an accurate simulation of the same normal form.
        fa = dyn_a.f
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
