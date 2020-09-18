module BDD

using CUDD

## Helper functions

_One(dd::Ptr{Manager}) = CUDD.Cudd_ReadOne(dd)
_Zero(dd::Ptr{Manager}) = CUDD.Cudd_ReadLogicZero(dd)
_Deref(dd::Ptr{Manager}, node::Ptr{Node}) = Cudd_RecursiveDeref(dd, node)
_Ref(node::Ptr{Node}) = CUDD.Cudd_Ref(node)
_IthVar(dd::Ptr{Manager}, idx::Cint) = CUDD.Cudd_bddIthVar(dd, idx)
# Remember that indexing starts at zero in C
_Size(dd::Ptr{Manager}) = Cint(CUDD.Cudd_ReadSize(dd))
_AND(dd::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddAnd(dd, f, g)
_OR(dd::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddOr(dd, f, g)
_XOR(dd::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddXor(dd, f, g)
_NOT(dd::Ptr{Manager}, f::Ptr{Node}) = _XOR(dd, f, _One(dd))
_Eval(dd::Ptr{Manager}, f::Ptr{Node}, values::Vector{Cint}) =
    CUDD.Cudd_Eval(dd, f, values) === _One(dd)

# Computes cube by providing indices of the variables instead of the variables
# themselves. The implementation is exactly the same as CUDD.Cudd_IndicesToCube
# except that we can also specify the phases.
function _Cube(dd::Ptr{Manager}, indices::Vector{Cint}, phases::Vector{Cint})
    @assert length(indices) === length(phases)
    n = length(indices)
    cube = _One(dd); _Ref(cube)
    for r in n:-1:1 # In CUDD it is implemented in backward order, so it mimiched...
        var = _IthVar(dd, indices[r])
        if phases[r] === Cint(1)
            tmp = _AND(dd, var, cube); _Ref(tmp)
        elseif phases[r] === Cint(0)
            not_var = _NOT(dd, var); _Ref(not_var)
            tmp = _AND(dd, not_var, cube); _Ref(tmp)
            _Deref(dd, not_var)
        else
            throw("Phases must be a Cint equal to 0 or 1")
        end
        _Deref(dd, cube)
        cube = tmp
    end
    _Deref(dd, cube)
    return cube
end

## BDD + vector to feed values without allocating too much
mutable struct BDDManager
    dd::Ptr{Manager}
    values::Vector{Cint}
end

function BDDManager()
    mng = BDDManager(CUDD.Cudd_Init(), Cint[])
    finalizer(mng -> CUDD.Cudd_Quit(mng.dd), mng)
    return mng
end

function add_var!(mng::BDDManager)
    dd = mng.dd
    nvars = _Size(dd)
    newvar = _IthVar(dd, nvars)
    resize!(mng.values, nvars + 1)
    return (newvar, nvars)
end

## A set a BDD variables depending on a manager.
# The goal is that BDD functions defined with a given Cluster as support are
# automatically updated when a variable is added to the support.
# The updating rule for the moment is:
# {new f}(x_1,...,x_n,x_{n+1}) = f(x_1,...,x_n) ∧ ¬x_{n+1}
# TODO: Add other updating rule depending on the type of the function (create new
# types for that), or either on the function itself (via an updating function
# stored in struct along with the boolean function).
mutable struct VariablesCluster
    mng::BDDManager
    nvars::Int
    indices::Vector{Cint}
    phases::Vector{Cint}
    Rroots::Vector{Ref{Ptr{Node}}}
end

function VariablesCluster(mng::BDDManager)
    return VariablesCluster(mng, 0, Cint[], Cint[], Ref{Ptr{Node}}[])
end

function add_var!(vc::VariablesCluster)
    dd = vc.mng.dd
    newvar, nvars = add_var!(vc.mng)
    resize!(vc.mng.values, nvars + 1)
    push!(vc.indices, nvars)
    vc.nvars += 1
    resize!(vc.phases, vc.nvars)
    for Rroot in vc.Rroots
        not_var = _NOT(dd, newvar); _Ref(not_var)
        tmp = _AND(dd, Rroot[], not_var); _Ref(tmp)
        _Deref(dd, not_var); _Deref(dd, Rroot[])
        Rroot[] = tmp
    end
end

function set_values!(vc)
    for (idx, phase) in zip(vc.indices, vc.phases)
        vc.mng.values[idx+1] = phase
    end
end

## -
include("integerset.jl")

end # module
