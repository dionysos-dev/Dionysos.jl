export EllipsoidalBackwardConfig,
    EllipsoidalBackwardRefinementOptions, BackwardStepRecord, EllipsoidalCertificationResult

struct EllipsoidalBackwardConfig{TX, TU, TW, Backend, Opts}
    Xdom::TX
    Udom::TU
    Wdom::TW
    backend::Backend
    options::Opts
end

Base.@kwdef struct EllipsoidalBackwardRefinementOptions
    state_scaling::Union{Nothing, Vector{Float64}} = nothing
end

struct BackwardStepRecord{TE, TK, TS}
    k::Int
    status::Symbol
    cost::Union{Nothing, Float64}
    ellipsoid::TE
    kappa::TK
    summary::TS
end

function BackwardStepRecord(
    k::Integer,
    status::Symbol,
    cost,
    ellipsoid,
    kappa,
    summary::TS,
) where {TS}
    return BackwardStepRecord(
        Int(k),
        status,
        cost === nothing ? nothing : Float64(cost),
        ellipsoid,
        kappa,
        summary,
    )
end

struct EllipsoidalCertificationResult{S, CTRL, LMI}
    success::Bool
    failed_k::Union{Nothing, Int}
    solve_time_sec::Float64
    steps::Vector{S}
    controller::CTRL
    lmi_data::LMI
end

function EllipsoidalCertificationResult(
    success::Bool,
    failed_k::Union{Nothing, Integer},
    solve_time_sec::Real,
    steps::Vector{S},
    controller,
    lmi_data,
) where {S}
    return EllipsoidalCertificationResult(
        success,
        failed_k === nothing ? nothing : Int(failed_k),
        Float64(solve_time_sec),
        steps,
        controller,
        lmi_data,
    )
end
