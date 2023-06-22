module TestMain
using Test

using Random
using LinearAlgebra
using StaticArrays
using IntervalArithmetic
using JuMP, SDPA

using Dionysos
const DI = Dionysos
const ST = DI.System
const UT = DI.Utils
const SY = DI.Symbolic

const TOL = 1e-6
const SEED = 0
Random.seed!(SEED)

@testset "ConstantController" begin
    c = [0.1; 0.2]
    cont = ST.ConstantController(c)
    @test cont.c === c
    @test (cont.c_eval)(1.) === c
end

@testset "check_feasibility" begin
    # Inspired by examples/SFRRT/single_transition.jl

    sm(M) = SMatrix{size(M, 1), size(M, 2)}(M)
    sv(M) = SVector{size(M, 1)}(M)

    Ac = sm([0.0 1.0; 1.0 -1.0])
    Bc = sm([1.0; 1.0])
    gc = sv(zeros(2, 1))
    E = sm([1.0; 1.0])

    nx = 2; nu = 1; nw = 1

    xbar = zeros(nx); ubar = zeros(nu); wbar = zeros(nw)
    ΔX = IntervalBox(-1.0 .. 1.0, nx)
    ΔU = IntervalBox(-10 * 2 .. 10 * 2, nu)
    ΔW = IntervalBox(-0.0 .. 0.0, nw)
    X̄ = IntervalBox(xbar .+ ΔX)
    Ū = IntervalBox(ubar .+ ΔU)
    W̄ = IntervalBox(wbar .+ ΔW)
    L = [0.0; 0.0; 0.0]

    Ubound = 10
    Ub = IntervalBox(-Ubound .. Ubound, nu)

    Uaux = diagm(1:nu)
    U = [(Uaux .== i) ./ Ub[i].hi for i in 1:nu]

    E1 = UT.Ellipsoid(Matrix{Float64}(I(nx)), [0.0; 0.0])
    E2 = UT.Ellipsoid(Matrix{Float64}(I(nx)) * (1 / 15), [3.0; 3.0])

    sys = ST.AffineApproximationDiscreteSystem(Ac, Bc, gc, E, X̄, Ū, W̄, L)

    W = 100. * [-1 1]
    S = Matrix{Float64}(I(nx + nu + 1))

    optimizer = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
    _, cont, _ =
        SY.transition_fixed(sys.constrainedAffineSys, E1, E2, U, W, S, optimizer)

    f_eval = ST.get_f_eval(sys)
    
    c_eval = ST.get_c_eval(cont)
    @test ST.check_feasibility(E1, E2, f_eval, c_eval, nw, Ub; N = 100)
end

end