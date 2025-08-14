module TestAlloc

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
using LinearAlgebra

# Test helper function to build example system
function build_test_system(; n_per_dim = 20, tstep = 1.0, input_step = 1.0)
    # 3D state domain
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)
    Xgrid = DO.GridFree(lb, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

    # 3D input domain
    lb_u = SVector(-1.0, -1.0, -1.0)
    ub_u = SVector(1.0, 1.0, 1.0)
    h_u = SVector(input_step, input_step, input_step)
    Ugrid = DO.GridFree(lb_u, h_u)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb_u, ub_u), DO.OUTER)

    # Linear dynamics: dx/dt = A*x + B*u
    A = @SMatrix [
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        -1.0 -1.0 -1.0
    ]
    B = @SMatrix [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0
    ]
    function F_sys(x, u)
        return A * x + B * u
    end

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        3,
        3,
        nothing,
        nothing,
    )

    return Xfull, Ufull, concrete_system
end

function _test_alloc(abstract_system, concrete_system_approx)
    inputs = @inferred Dionysos.Symbolic.enum_inputs(abstract_system)
    abstract_input = @inferred first(inputs)
    # FIXME `symmodel.metadata[:original_symmodel]` is type-unstable since `metadata` value type is `Any`
    #concrete_input = @infered Dionysos.Symbolic.get_concrete_input(abstract_system, abstract_input)
    concrete_input = Dionysos.Symbolic.get_concrete_input(
        abstract_system,
        abstract_input,
    )::SVector{3, Float64}
    abstract_state = @inferred first(Dionysos.Symbolic.enum_states(abstract_system))
    concrete_state =
        @inferred Dionysos.Symbolic.get_concrete_state(abstract_system, abstract_state)
    system_map = @inferred ST.get_system_map(concrete_system_approx)
    # Compilation
    system_map(concrete_state, concrete_input)
    @test 0 == @allocated system_map(concrete_state, concrete_input)
end

function test_alloc()
    Xfull, Ufull, concrete_system = build_test_system(; n_per_dim = 2)

    # Create discrete system
    cont_center = ST.ContinuousTimeCenteredSimulation(concrete_system)
    discrete_system = ST.discretize(cont_center, 1.0)

    sym_serial = SY.NewSymbolicModelListList(Xfull, Ufull)
    return _test_alloc(sym_serial, discrete_system)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestAlloc.runtests()
