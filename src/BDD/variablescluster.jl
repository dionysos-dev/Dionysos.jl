mutable struct VariablesCluster
    mng::BDDManager
    indices_::Vector{Cint}
    auxindices_::Vector{Cint}
    phases1_::Vector{Cint}
    phases2_::Vector{Cint}
    auxphases_::Vector{Cint}
    roots_::Vector{Ptr{Node}}
end

function add_var!(varcl::VariablesCluster)
    nvars = _Size(mng.dd)
    _IthVar(mng, nvars)
    push!(indices, nvars)
end

function add_var!(varcl::VariablesCluster)
    nvars = _Size(mng.dd)
    _IthVar(mng, nvars)
    push!(indices, nvars)
    push!(auxindices, nvars)
end
