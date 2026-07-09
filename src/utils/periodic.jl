# Periodic-domain coordinate wrapping (single home; used by System's closed-loop
# simulation and Mapping's periodic grid mapping).

"Wrap the scalar `x` into `[start, start + period)`."
wrap_value(x, start, period) = mod(x - start, period) + start

"""
    wrap_coord(x::SVector{N, T}, periodic_dims::SVector{P, Int}, periods::SVector{P, T}; start = zeros(SVector{P, T}))

Wraps the vector `x` into a periodic domain along dimensions specified in `periodic_dims`,
with period lengths `periods` and optional offset `start`.

# Arguments
- `x`: The coordinate vector to wrap.
- `periodic_dims`: Indices of the periodic dimensions.
- `periods`: Period lengths for the periodic dimensions.
- `start` (optional): Starting values of the periodic domains (defaults to `0.0`).

# Returns
A wrapped `SVector` where each periodic dimension is mapped to the interval `[start[i], start[i] + periods[i])`.
"""
function wrap_coord(
    x::SVector{N, T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {N, P, T}
    return SVector{N, T}(
        ntuple(d -> begin
            i = findfirst(isequal(d), periodic_dims)
            i === nothing ? x[d] : wrap_value(x[d], start[i], periods[i])
        end, N),
    )
end

function get_periodic_wrapper(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {P, T}
    return x -> wrap_coord(x, periodic_dims, periods; start = start)
end
