module Allocation

using Test
using StaticArrays
using LinearAlgebra
using MathematicalSystems
using Dionysos

const DI = Dionysos
const MP = DI.Mapping
const ST = DI.System
const SY = DI.Symbolic

# -----------------------
# Helper: build finite explicit grid mappings (3D)
# -----------------------
function build_test_mappings(; n_per_dim::Int = 2, input_step::Float64 = 1.0)
    # ---- X mapping: positions in a small 3D box ----
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)

    Xgrid = MP.GridFree(lb, h)
    Xmap = MP.ExplicitGridMapping(Xgrid)

    # enumerate positions 0:(n_per_dim-1) in each dimension
    for i in 0:(n_per_dim - 1), j in 0:(n_per_dim - 1), k in 0:(n_per_dim - 1)
        MP.add_pos!(Xmap, (i, j, k))
    end

    # ---- U mapping: small box [-1,1]^3 sampled with step input_step ----
    lb_u = SVector(-1.0, -1.0, -1.0)
    h_u = SVector(input_step, input_step, input_step)
    Ugrid = MP.GridFree(lb_u, h_u)
    Umap = MP.ExplicitGridMapping(Ugrid)

    # positions corresponding to coords -1, 0, 1 when step=1
    # (pos indices here depend on your GridFree convention; we just insert a few explicitly)
    for i in 0:2, j in 0:2, k in 0:2
        MP.add_pos!(Umap, (i, j, k))
    end

    return Xmap, Umap
end

# -----------------------
# Helper: build example continuous system
# -----------------------
function build_test_system()
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

    F_sys(x, u) = A * x + B * u

    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        3,
        3,
        nothing,
        nothing,
    )
end

function _test_alloc(symmodel, discrete_system)
    inputs = @inferred SY.enum_inputs(symmodel)
    abstract_input = @inferred first(inputs)

    # Make the type stable explicitly (older code comment still applies sometimes)
    concrete_input = SY.get_concrete_input(symmodel, abstract_input)::SVector{3, Float64}

    abstract_state = @inferred first(SY.enum_states(symmodel))
    concrete_state = @inferred SY.get_concrete_state(symmodel, abstract_state)

    system_map = @inferred ST.get_system_map(discrete_system)

    # compile once
    system_map(concrete_state, concrete_input)

    # allocation check
    @test 0 == @allocated system_map(concrete_state, concrete_input)
end

function test_alloc()
    Xmap, Umap = build_test_mappings(; n_per_dim = 2, input_step = 1.0)
    concrete_system = build_test_system()

    cont_center = ST.ContinuousTimeCenteredSimulation(concrete_system)
    discrete_system = ST.discretize(cont_center, 1.0)

    # By default SymbolicModelList uses MappingSet{N}() which is "all states of mapping"
    sym = SY.SymbolicModelList(Xmap, Umap)
    return _test_alloc(sym, discrete_system)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith(String(name), "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end # module

Allocation.runtests()
