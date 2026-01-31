module Artifacts
using JLD2

export Artifact, save_artifact, load_artifact

struct Artifact
    kind::Symbol
    version::Int
    meta::Dict{String, Any}
    payload::Dict{String, Any}
end

function save_artifact(path::AbstractString, art::Artifact)
    @save path art
    return path
end

function load_artifact(path::AbstractString)::Artifact
    art = nothing
    @load path art
    return art::Artifact
end

end # module
