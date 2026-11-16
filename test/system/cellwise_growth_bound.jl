module TestMain

import Dionysos
import LazySets
import MathOptInterface as MOI
import Random
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))
import MathematicalSystems as MS
const UGA = AB.UniformGridAbstraction

# `ST.CellwiseBound`: a growth bound whose Jacobian bound is taken over an enclosure of the step
# from each cell. Test system: the pendulum ẋ₁ = x₂, ẋ₂ = −sin x₁ + u, whose Jacobian
# [0 1; −cos x₁ 0] has the Metzler bound [0 1; max|cos x₁| 0] over a box.
pendulum(x, u) = SVector(x[2], -sin(x[1]) + u[1])

lo(cell, i) = LazySets.low(cell, i)
hi(cell, i) = LazySets.high(cell, i)

# A crude but valid enclosure of the step: x₂ moves by at most t (1 + |u|), and x₁ by at most t
# times the largest |x₂| on the way.
function enclosure(cell, u, τ)
    w = τ * (1 + abs(u[1]))
    E2 = (lo(cell, 2) - w, hi(cell, 2) + w)
    m2 = max(abs(E2[1]), abs(E2[2]))
    return ((lo(cell, 1) - τ * m2, hi(cell, 1) + τ * m2), E2)
end

function max_abs_cos(a)
    k = ceil(a[1] / π)
    k * π <= a[2] && return 1.0 # a multiple of π, where |cos| = 1, lies in the interval
    return max(abs(cos(a[1])), abs(cos(a[2])))
end
bound(E, _) = SMatrix{2, 2}(0.0, max_abs_cos(E[1]), 1.0, 0.0)
global_bound(_) = SMatrix{2, 2}(0.0, 1.0, 1.0, 0.0)

system = MS.ConstrainedBlackBoxControlContinuousSystem(
    pendulum,
    2,
    1,
    LazySets.Hyperrectangle(; low = SVector(-π, -2.0), high = SVector(π, 2.0)),
    LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0)),
)
const τ = 0.1

discrete(jacobian_bound) = ST.discretize(
    ST.ContinuousTimeGrowthBound(system; jacobian_bound, ngrowthbound = 10),
    τ;
    num_substeps = 10,
)
reach(approx, cell, u) =
    ST.reach_set(approx, cell, u, ST.input_cache(approx, cell.radius, u))

@testset "a constant cellwise bound is the per-input bound" begin
    cellwise = discrete(ST.CellwiseBound((cell, u, t) -> cell, (E, u) -> global_bound(u)))
    per_input = discrete(global_bound)
    cell = LazySets.Hyperrectangle(SVector(0.3, -0.4), SVector(0.05, 0.05))
    for u in (SVector(-1.0), SVector(0.0), SVector(1.0))
        @test reach(cellwise, cell, u) == reach(per_input, cell, u)
    end
end

@testset "sound and tighter than the global bound" begin
    cellwise = discrete(ST.CellwiseBound(enclosure, bound))
    per_input = discrete(global_bound)
    rng = Random.MersenneTwister(1)
    fine(x, u) = ST.runge_kutta4(pendulum, x, u, τ, 200)
    for _ in 1:50
        c = SVector(2π * rand(rng) - π, 3 * rand(rng) - 1.5)
        cell = LazySets.Hyperrectangle(c, SVector(0.05, 0.05))
        u = SVector(rand(rng, (-1.0, 0.0, 1.0)))
        box = reach(cellwise, cell, u)
        @test all(box.radius .<= reach(per_input, cell, u).radius)
        for _ in 1:20
            x = c .+ (2 .* SVector(rand(rng), rand(rng)) .- 1) .* cell.radius
            @test fine(x, u) in box
        end
    end
end

@testset "no enclosure disables the input" begin
    picky =
        ST.CellwiseBound((cell, u, t) -> u[1] > 0 ? nothing : enclosure(cell, u, t), bound)
    @test reach(
        discrete(picky),
        LazySets.Hyperrectangle(SVector(0.0, 0.0), SVector(0.05, 0.05)),
        SVector(1.0),
    ) === nothing

    function abstraction(jacobian_bound)
        X = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
        sys = MS.ConstrainedBlackBoxControlContinuousSystem(pendulum, 2, 1, X, system.U)
        target =
            LazySets.Hyperrectangle(; low = SVector(-0.1, -0.1), high = SVector(0.1, 0.1))
        start = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(0.6, 0.6))
        problem = PR.OptimalControlProblem(
            sys,
            start,
            target,
            nothing,
            (x, u) -> 1.0,
            PR.Infinity(),
        )
        opt = MOI.instantiate(UGA.Optimizer)
        set(k, v) = MOI.set(opt, MOI.RawOptimizerAttribute(k), v)
        set("concrete_problem", problem)
        set("state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1)))
        set("input_grid", MP.GridFree(SVector(0.0), SVector(1.0)))
        set("time_step", τ)
        set("approx_mode", UGA.GROWTH)
        set("jacobian_bound", jacobian_bound)
        set("print_level", 0)
        MOI.optimize!(opt)
        abs_sys = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
        input(s) = SY.get_concrete_input(abs_sys, s)[1]
        succ = Dict{Tuple{Int, Float64}, Set{Int}}()
        for (t, q, s) in SY.enum_transitions(SY.get_automaton(abs_sys))
            push!(get!(succ, (q, input(s)), Set{Int}()), t)
        end
        return succ
    end
    cellwise = abstraction(picky)
    per_input = abstraction(global_bound)
    @test !isempty(cellwise)
    @test all(u <= 0 for (_, u) in keys(cellwise))
    @test any(u > 0 for (_, u) in keys(per_input))
    @test all(haskey(per_input, k) && s ⊆ per_input[k] for (k, s) in cellwise)
end

end # module
