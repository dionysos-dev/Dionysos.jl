module TestControllers

using Test
import MathematicalSystems as MS

using Dionysos
const DI = Dionysos
const ST = DI.System

# Minimal protocol implementations to exercise the trait defaults.
struct ToyStatic <: ST.AbstractController end
ST.controller_kind(::ToyStatic) = ST.StaticKind()
ST.output_control(::ToyStatic, mem, y) = 2 * y

struct ToyDynamic <: ST.AbstractController end
ST.controller_kind(::ToyDynamic) = ST.DynamicKind()

@testset "Trait defaults" begin
    c = ToyStatic()
    @test ST.initial_state(c) === nothing
    @test ST.update_state(c, nothing, 1.0) === nothing
    @test ST.domain(c) === nothing
    @test ST.is_defined(c, nothing, 1.0)
    @test ST.output_control(c, nothing, 3.0) == 6.0

    d = ToyDynamic()
    @test_throws ErrorException ST.initial_state(d)
    @test_throws ErrorException ST.update_state(d, nothing, 1.0)
end

@testset "ControlTable + DiscreteStaticController" begin
    tab = ST.ControlTable(3)
    ST.add_control!(tab, 1, 4)
    ST.add_control!(tab, 1, 5)
    ST.set_control!(tab, 2, 7)
    @test tab(1) == [4, 5]
    @test tab(2) == [7]
    @test isempty(tab(3))

    ctrl = ST.DiscreteStaticController(Set([1, 2]), tab, false)
    @test ST.controller_kind(ctrl) == ST.StaticKind()
    @test ST.initial_state(ctrl) === nothing
    @test ST.is_defined(ctrl, nothing, 1)
    @test ST.output_control(ctrl, nothing, 1) == 4
    @test ST.output_control(ctrl, nothing, 2) == 7
    @test !ST.is_defined(ctrl, nothing, 3)     # empty input list
    @test ST.output_control(ctrl, nothing, 3) === nothing
end

@testset "AffineController / as_controller" begin
    K = [1.0 0.0; 0.0 2.0]
    b = [0.5, -0.5]
    result = ST.TransitionResult(true, MS.AffineMap(K, b), 1.0, nothing)
    ctrl = ST.as_controller(result)
    @test ctrl isa ST.AffineController
    @test ST.controller_kind(ctrl) == ST.StaticKind()
    x = [1.0, 1.0]
    @test ST.output_control(ctrl, nothing, x) == K * x + b

    @test ST.as_controller(ST.TransitionResult(false, nothing, nothing, nothing)) ===
          nothing
end

@testset "AutomatonMemoryController protocol" begin
    # Two abstract states; state 1 labeled :a. Memory: 1 --:a--> 2 (accepting),
    # stays otherwise. Product states: (qs, qa) -> p; product controller plays
    # input 9 from p = 1 and input 8 from p = 2.
    tab = ST.ControlTable(2)
    ST.add_control!(tab, 1, 9)
    ST.add_control!(tab, 2, 8)
    product_ctrl = ST.DiscreteStaticController(Set([1, 2]), tab, false)

    label_of_state = [2, 1]                      # qs = 1 -> label 2 (:a), qs = 2 -> empty
    step_map = Dict((1, 1) => 1, (1, 2) => 2, (2, 1) => 2, (2, 2) => 2)
    pid = Dict((1, 1) => 1, (2, 2) => 2)

    ctrl = ST.AutomatonMemoryController(1, label_of_state, 1, step_map, pid, product_ctrl)

    @test ST.controller_kind(ctrl) == ST.DynamicKind()
    qa = ST.initial_state(ctrl)
    @test qa == 1

    @test ST.is_defined(ctrl, qa, 1)
    @test ST.output_control(ctrl, qa, 1) == 9
    qa2 = ST.update_state(ctrl, qa, 1)           # sees :a -> memory 2
    @test qa2 == 2

    @test ST.output_control(ctrl, qa2, 2) == 8
    @test ST.update_state(ctrl, qa2, 2) == 2

    # (qs, qa) pair outside the product: undefined
    @test !ST.is_defined(ctrl, qa, 2)
    @test ST.output_control(ctrl, qa, 2) === nothing

    # measurement outside the labeled range falls back to the default label
    @test ST.update_state(ctrl, 1, 5) == 1
end

end # module
