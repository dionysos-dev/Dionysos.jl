module TestMain

using Test, LinearAlgebra
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const OP = DI.Optim
const AB = OP.Abstraction

sleep(0.1) # used for good printing
println("Started test")

function test_symbolic_controller(ctrl::ST.SymbolicController)
    ST.add_control!(ctrl, 1, 5)
    ST.add_control!(ctrl, 1, 6)

    @test ST.is_defined(ctrl, 1)
    syms = ST.get_all_controls(ctrl, 1)
    @test length(syms) == 2
    @test 5 ∈ syms && 6 ∈ syms
    @test ST.get_control(ctrl, 1) ∈ syms

    @test !ST.is_defined(ctrl, 2)
    @test 2 ∉ ST.domain(ctrl)
end

@testset "SymbolicControllerList" begin
    ctrl = ST.SymbolicControllerList()
    test_symbolic_controller(ctrl)
end
@testset "SymbolicControllerDict" begin
    ctrl = ST.SymbolicControllerDict()
    test_symbolic_controller(ctrl)
end

@testset "ConstantController" begin
    u = [1.0, 2.0]
    ctrl = ST.ConstantController(u)
    x = [0.0, 0.0]
    @test ST.get_control(ctrl, x) == u
    @test ST.is_defined(ctrl, x)
    @test ST.get_all_controls(ctrl, x) == [u]
end

@testset "AffineController" begin
    K = [1.0 0.0; 0.0 1.0]
    c = [0.0, 0.0]
    ℓ = [1.0, 1.0]
    ctrl = ST.AffineController(K, c, ℓ)
    x = [2.0, 3.0]
    @test ST.get_control(ctrl, x) ≈ x + ℓ
    @test ST.is_defined(ctrl, x)
end

@testset "BlackBoxContinuousController" begin
    f(x) = x .* 2
    isdef(x) = norm(x) < 10
    ctrl = ST.BlackBoxContinuousController(f, isdef)

    x1 = [1.0, 2.0]
    x2 = [100.0, 200.0]

    @test ST.is_defined(ctrl, x1)
    @test !ST.is_defined(ctrl, x2)
    @test ST.get_control(ctrl, x1) == f(x1)
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
