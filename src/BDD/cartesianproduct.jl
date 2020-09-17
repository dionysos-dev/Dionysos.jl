mutable struct CartesianProduct
    mng::Ptr{Manager}
    nslices::Int
    indices_::Vector{Vector{Cint}}
    auxindices_::Vector{Vector{Cint}}
    phases1_::Vector{Vector{Cint}}
    phases2_::Vector{Vector{Cint}}
    auxphases_::Vector{Vector{Cint}}
    values::Vector{Cint}
    roots_::Vector{Vector{Ptr{Node}}}
end

function CartesianProduct(nslices::Int)
    mng = CUDD.Cudd_Init()
    ARGS = ntuple(k -> [Cint[] for i in 1:nslices], 5)
    roots_ = [Ptr{Node}[] for i in 1:nslices]
    cp = CartesianProduct(mng, ARGS..., Cint[], roots_)
    finalizer(cp -> CUDD.Cudd_Quit(cp.mng), cp)
    return cp
end

function add_slice!(cp::CartesianProduct)
    append!(cp.indices_, Cint[])
    append!(cp.auxindices_, Cint[])
    append!(cp.phases1_, Cint[])
    append!(cp.phases2_, Cint[])
    append!(cp.auxphases_, Cint[])
    append!(cp.roots_, Ptr{Node}[])
    cp.nslices += 1
    return cp
end

function add_root(cp::CartesianProduct, root, slices)
    for i in slices
        push!(cp.roots_[i], root)
    end
    return cp
end
