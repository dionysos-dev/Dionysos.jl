module RuntimeAPI

export AbstractRuntimeController, step!, reset!

abstract type AbstractRuntimeController end

function step!(ctrl::AbstractRuntimeController; x, t=0.0, kwargs...)
    error("step! not implemented for $(typeof(ctrl))")
end

reset!(ctrl::AbstractRuntimeController; kwargs...) = (; ok=true)

end # module
