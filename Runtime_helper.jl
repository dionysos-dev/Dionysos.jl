Base.@kwdef mutable struct PlannerBackend
    tokenizer
    model
    device
    prefix::String
    max_input_length::Int
    max_target_length::Int
    abstract_system
    ts::TSFromAdj
    Label
    inner_sets::Dict{Symbol, Set{Int}}
    x0_cont::SVector{3,Float64}
    x0_abs_start::Int
    hx::SVector{3,Float64}
end

Base.@kwdef mutable struct BuchiBuildResult
    nl::String
    formula::String
    aut::Automaton
    labels_ts::Vector{Set{Symbol}}
end

