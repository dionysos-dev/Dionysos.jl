# Expose ellipsoidal-backward public data types.
export EllipsoidalBackwardConfig, BackwardStepRecord, EllipsoidalCertificationResult

# Store static configuration for ellipsoidal backward certification.
struct EllipsoidalBackwardConfig{TX, TU, TW, Backend, Opts}
    # State domain descriptor.
    Xdom::TX
    # Input domain descriptor.
    Udom::TU
    # Disturbance domain descriptor.
    Wdom::TW
    # Backend descriptor for optimization/LMI execution.
    backend::Backend
    # Additional user options container.
    options::Opts
end

# Store one backward-step output record.
struct BackwardStepRecord{TE, TK, TS}
    # Backward index k.
    k::Int
    # Step status tag.
    status::Symbol
    # Optional step cost.
    cost::Union{Nothing, Float64}
    # Optional ellipsoid object.
    ellipsoid::Union{Nothing, TE}
    # Optional local control law object.
    kappa::Union{Nothing, TK}
    # Extensible summary payload.
    summary::TS
end

# Store global certification output.
struct EllipsoidalCertificationResult{S, CTRL, LMI}
    # Global success flag.
    success::Bool
    # First failing backward index, if any.
    failed_k::Union{Nothing, Int}
    # Certification solve time in seconds.
    solve_time_sec::Float64
    # Ordered list of backward-step records.
    steps::Vector{S}
    # Optional synthesized controller payload.
    controller::Union{Nothing, CTRL}
    # Optional raw LMI payload.
    lmi_data::Union{Nothing, LMI}
end

function EllipsoidalCertificationResult(
    success::Bool,
    failed_k::Union{Nothing, Int},
    solve_time_sec::Float64,
    steps::Vector{S},
    controller::Nothing,
    lmi_data::Nothing,
) where {S}
    return EllipsoidalCertificationResult{S, Nothing, Nothing}(
        success,
        failed_k,
        solve_time_sec,
        steps,
        controller,
        lmi_data,
    )
end

function EllipsoidalCertificationResult(
    success::Bool,
    failed_k::Union{Nothing, Int},
    solve_time_sec::Float64,
    steps::Vector{S},
    controller::Nothing,
    lmi_data::LMI,
) where {S, LMI}
    return EllipsoidalCertificationResult{S, Nothing, LMI}(
        success,
        failed_k,
        solve_time_sec,
        steps,
        controller,
        lmi_data,
    )
end

function EllipsoidalCertificationResult(
    success::Bool,
    failed_k::Union{Nothing, Int},
    solve_time_sec::Float64,
    steps::Vector{S},
    controller::CTRL,
    lmi_data::Nothing,
) where {S, CTRL}
    return EllipsoidalCertificationResult{S, CTRL, Nothing}(
        success,
        failed_k,
        solve_time_sec,
        steps,
        controller,
        lmi_data,
    )
end
