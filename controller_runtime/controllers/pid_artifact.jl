module PIDArtifact

using StaticArrays
import ..Artifacts: Artifact

export make_pid_artifact

function make_pid_artifact(;
    tstep::Float64,
    Kp,
    Ki,
    Kd,
    ref,
    umin = nothing,
    umax = nothing,
    antiwindup::Bool = true,
    meta::Dict{String,Any} = Dict{String,Any}(),
)
    payload = Dict{String,Any}(
        "tstep" => tstep,
        "Kp"    => Kp,
        "Ki"    => Ki,
        "Kd"    => Kd,
        "ref"   => ref,
        "umin"  => umin,
        "umax"  => umax,
        "antiwindup" => antiwindup,
        "e0"    => SVector(0.0, 0.0),
        "I0"    => SVector(0.0, 0.0),
    )

    meta2 = copy(meta)
    meta2["controller_kind"] = "pid"

    return Artifact(:pid, 1, meta2, payload)
end

end # module
