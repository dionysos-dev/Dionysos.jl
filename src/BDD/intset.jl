"""
    mutable struct IntSet <: AbstractSet{<:Integer}

Same as `Base.Set{<:Integer}` but with `CUDD`.
"""
mutable struct IntSet{T<:Integer} <: AbstractSet{T}
    mng::Ptr{Manager}
    root::Ptr{Node}
    indices_::Tuple{Vector{Cint}}
    auxindices_::Tuple{Vector{Cint}}
    phases1_::Tuple{Vector{Cint}}
    phases2_::Tuple{Vector{Cint}}
    auxphases_::Tuple{Vector{Cint}}
    values::Vector{Cint}
end

function IntSet{T}() where T
    mng = CUDD.Cudd_Init()
    ARGS = ntuple(k -> tuple(Cint[]), 5)
    root = _Zero(mng); _Ref(root)
    set = IntSet{T}(mng, root, ARGS..., Cint[])
    finalizer(set) do set
        CUDD.Cudd_Quit(set.mng)
    end
    return set
end
IntSet() = IntSet{Int}()

Base.eltype(::Type{IntSet{T}}) where T = T
Base.empty(::IntSet, ::Type{T}=Int) where T = IntSet{T}()
Base.emptymutable(::IntSet, ::Type{T}=Int) where T = IntSet{T}()

function _phases1!(set::IntSet, x)
    indices = set.indices_[1]
    auxindices = empty!(set.auxindices_[1])
    phases = set.phases1_[1]
    auxphases = empty!(set.auxphases_[1])
    _compute_phases!(set.mng, phases, indices, auxphases, auxindices, x)
end

#=
function _phases2!(set::IntSet, x, y)
    phases1 = set.phases1
    phases2 = set.phases2
    newphases1 = empty!(set.newphases1)
    newphases2 = empty!(set.newphases2)
    newvars = empty!(set.newvars)
    auxphases = empty!(set.auxphases)
    nvars = length(set.variables)
    indices = 1:nvars
    _compute_phases!(set.mng, phases1, newvars, newphases1, indices, x)
    nnewvars1 = length(newphases1)
    resize!(phases2, nvars + nnewvars1)
    indices = 1:(nvars + nnewvars1)
    _compute_phases!(set.mng, phases2, newvars, newphases2, indices, y)
    nnewvars2 = length(newphases2)
    for i = 1:nnewvars1
        push!(phases1, newphases1[i])
    end
    for i = 1:nnewvars2
        push!(phases1, Cint(0))
        push!(phases2, newphases2[i])
    end
    append!(set.variables, newvars)
    resize!(auxphases, nnewvars1 + nnewvars2); fill!(auxphases, Cint(0))
end

# Returns the smallest `rank` such that `phases1[indices[r]] === phases2[indices[r]]`
# for all `r` ≧ `rank`
function _common_part(phases1, phases2, indices)
    nidx = length(indices)
    for r in nidx:-1:1
        phases1[indices[r]] !== phases2[indices[r]] && return r
    end
    return 0
end

function _upper_bound(root, phases, indices, variables, r)

end


function push_interval!(set::IntSet{UInt}, x::UInt, y::UInt)
    x > y && return set
    _phases2!(set::IntSet, x, y)
    c = _Cube(set.mng, set.newvars, set.auxphases); _Ref(c)
    tmp1 = CUDD.Cudd_bddAnd(set.mng, set.root, c); _Ref(tmp1)
    _Deref(set.mng, set.root); _Deref(set.mng, c)
    phases1 = set.phases1
    phases2 = set.phases2
    indices = 1:length(set.variables)
    rank = _common_part(phases1, phases2, indices)


    # c = _Cube(set.mng, set.variables, set.phases1); _Ref(c)
    # tmp2 = CUDD.Cudd_bddOr(set.mng, tmp1, c); _Ref(tmp2)
    # _Deref(set.mng, tmp1); _Deref(set.mng, c)
    # set.root = tmp2
    return set
end

# function _phase2(set::IntSet, x1, x2)
#     empty!(set.vars)
#     empty!(set.z_)
#     for idx in eachindex(set.variables)
#         set.phases_[idx] = _bit(x)
#         x >>>= 1
#     end
#     while x > 0
#         push!(set.phases_, _bit(x))
#         # As the `mng` is only used by this struct, `bddIthVar` should be
#         # the same as `bddNewVar`.
#         newvar = CUDD.Cudd_bddIthVar(set.mng, length(set.variables))
#         push!(set.variables, newvar)
#         push!(set.vars_, newvar)
#         push!(set.z_, zero(Cint))
#         x >>>= 1
#     end
# end
=#

# Return 0 if there are enough bits to represent x.
# Returning a bool would have been more straighforward but this is to be consistent
# with IntTupleSet.
function _phases1_trunc!(set::IntSet, x)
    not_trunc = _compute_phases_trunc!(set.phases1_[1], x)
    return not_trunc ? 0 : 1
end

function Base.iterate(set::IntSet{T}, state::T=zero(T)) where T
    _phases1_trunc!(set, state) > 0 && return nothing
    _eval_phases(set) && return (state, state + 1)
    return iterate(set, state + 1)
end
